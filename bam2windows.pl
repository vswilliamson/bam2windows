# Prepare input file for CNAnorm discrete normalization from a pair of sam/bam
# files. If you use bam format, samtools must be in your $PATH
# Type:
# perl <nameOfThisScript>
# to see usage
#
#
# Copyright (C) 2010 Stefano Berri <s.berri@leeds.ac.uk> or <s.berri@gmail.com>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Version 0.3
# created 25/08/2010
# last update 21/10/2011
my $version = "0.3.4";



use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;

my %par = setDefaultPars();
my $usage = setUsage(%par);

GetOptions (\%par, 'window=i', 'gc_file=s', 'readNum=i', 'genomeSize=i', 
    'qualityThreshold=i', 'tmpDir=s', 'testTemp', 'controlTemp', 
    'makeTempOnly', 'saveTest=s', 'saveControl=s', 'chrFile=s');

parameterCheck(%par);

if (@ARGV == 0 || @ARGV > 2){print $usage; exit}
if (@ARGV == 1 && !defined($par{'makeTempOnly'})){print $usage; exit}


my $testFile = shift;
my $controlFile = shift;
# load default chromosome length
# my %chrLength = chrLenght($par{'chrFile'});
my %chrLength = chrLength($testFile, $controlFile);
# create filter object
my $filters_hr = setFilters(%par);

my ($countTestFile, $readNumTest) = sam2simple($testFile, $par{'tmpDir'}, 
    $par{'testTemp'}, $filters_hr);
# exit if (!defined($controlFile) && $par{'makeTempOnly'});

my ($countControlFile, $readNumControl) = sam2simple($controlFile, $par{'tmpDir'},
    $par{'controlTemp'}, $filters_hr) if ($controlFile);
# exit if ($par{'makeTempOnly'});

# if we don't need only the temp files...
if (!defined($par{'makeTempOnly'})){
    # how many reads in the smallest of the two files?
    my $NumReads = min($readNumTest, $readNumControl);
    if (!$par{'window'}) {
        # set window size accordingly
        $par{'window'} = getWindowSize($NumReads, $par{'genomeSize'},
            $par{'readNum'});
    }


    # produce counts
    warn "Counting reads over window $par{'window'} bp wide...\n";
    my $hr = windowCount($par{'window'}, $countTestFile, $countControlFile, 
        $filters_hr, \%chrLength);
    
    if ($par{'gc_file'} ne '') {
        warn "  Calculating GC content...\n";
        addGC_content($hr, $par{'gc_file'}, $par{'window'});        
    }
    printHash($hr); 

}


# rename or delete temporary file
if (! renameOrDelete($countTestFile, $par{'saveTest'}, $par{'testTemp'})) {
    warn "Impossible to rename or delete test temporary file\n";
}

if ($controlFile){
    if (! renameOrDelete($countControlFile, $par{'saveControl'}, $par{'controlTemp'})) {
        warn "Impossible to rename or delete control temporary file\n";
    }
}





#### SUBS ####
sub parameterCheck {
    my %par = @_;
    if ($par{'makeTempOnly'}){
        if ((!$par{'saveTest'}) && (!$par{'saveControl'}) ){
            die "If --makeTempOnly is present --saveTest and/or --saveControl must be set too";
        }
    }
    if ($par{'chrFile'}){
        warn "Option 'chrFile' is deprecated and not used. sam/bam header is read instead\n";
    }
    return 0;
}


sub renameOrDelete {
    my ($file, $saveFile, $tempFile) = @_;
    my $success;
    if ($tempFile){
        $success = 1;
    } elsif ($saveFile ne ''){
        $success = rename $file, $saveFile;
    } else {
        $success = unlink ($file);
    }
    return $success;
}

sub min {
    my @nums = @_;
    my @sorted = sort {$a <=> $b} @nums;
    return $sorted[0];
}

sub sam2simple {
    my ($file, $dir, $tmpF, $filters_hr) = @_;
    my $numOfGoodReads = 0;
    
    if ($tmpF){
        open (TMP, $file) || die "Impossible to open temporary file $file\n";
        my $firstLine = <TMP>;
        $firstLine =~ /^# 0*(\d+)/;
        seek (TMP, $1, 0);
        
        my $found = 0;
        while (defined(my $line = <TMP>) && (! $found)){
            if ($line =~ /^# numOfGoodReads = (\d+)/){
                $numOfGoodReads = $1;
                $found = 1;
            }
        }
        if ($numOfGoodReads == 0){
            warn "Unable to find numOfGoodReads in the temporary file, calculating...\n";
            seek (TMP, 0, 0);
            while (<TMP>){
                if (! /^\s*#/){
                    $numOfGoodReads++;
                }
            }
        }
        return ($file, $numOfGoodReads);
    } else {
        
        my ($fileName) = fileparse($file);
        my $fh = getFileHandle($file);
        $dir =~ s/\/$//;
        my $tmpFile = $dir . "/" . createTmpFileName($fileName);
        # warn "creating temporary file $tmpFile...\n";
        open (OUT, ">$tmpFile") || die "Impossible to create temporary file $tmpFile: $!\n";
        # prepare space for pointer. Will be substitute at the end
        print OUT "# 000000000000\n";
        while (<$fh>){
            next if /^\@/;
            my @F = split(/\t/);
            my ($chr, $pos) = line2pos(\@F, $filters_hr);
            if (defined ($chr)) {
                print OUT "$chr\t$pos\n";
                $numOfGoodReads++;
            }
        }
        my $pos = sprintf("%012d", tell(OUT));
        print OUT "# numOfGoodReads = $numOfGoodReads\n";
        seek (OUT, 0, 0);
        print OUT "# $pos\n";
        
        close OUT;
        return ($tmpFile, $numOfGoodReads);
     }
}

sub createTmpFileName {
    my ($fileName) = @_;
    my $tmpFileName = $fileName ."_" . sprintf("%06x", int(rand(16777215))) . ".tmp";
    if ( -e $tmpFileName ) {
        print "recurring....\n";
        $tmpFileName = createTmpFileName($fileName);

    }
    return $tmpFileName;
}

sub getWindowSize {
    my ($TotNumReads, $genomeSize, $readPerWind)  = @_;
    my $windSize = $readPerWind * $genomeSize / $TotNumReads;
    return int($windSize);
}

sub getNumOfReads {
    my ($reads_vr, $testFile, $controlFile) = @_;
    if (${$reads_vr}) {
        return 1;
    } else {
        my $testR = 0;
        open (IN,  $testFile) || die "Impossible to open $testFile: $!\n";
        while(<IN>) {
            $testR++;
        }
        my $contrR = 0;
        open (IN,  $controlFile) || die "Impossible to open $controlFile: $!\n";
        while(<IN>) {
            $contrR++;
        }
        if ($testR < $contrR) {
            ${$reads_vr} = $testR;
        } else {
            ${$reads_vr} = $contrR;
        }
        close(IN);

    }
    return 1;
}

sub addGC_content {
    my ($hr, $inputFile, $wSize) = @_;
    
    if ($inputFile =~ /\.gz$/) {
        open (IN, "gzip -cd $inputFile |") || die "Impossible to open $inputFile: $!\n";
    } else {
        open (IN, $inputFile) || die "Impossible to open $inputFile: $!\n";
    }

    my $count_hr = {};
    my $chr;
    while(<IN>) {
        # in case it is windows or mac format
        s/\r//;
        chomp;
        # if (/chrom=chr(\S+)/) {
        if (/chrom=(\S+)/) {
            $chr = $1;
            $count_hr->{$chr} = {};
        } else {
            my ($pos, $GC) = split(/\t/, $_);
            my $start = int($pos/$wSize) * $wSize + 1;
#             if (ref ($count_hr->{$chr}->{$start}) ne 'ARRAY'){
#                 # initialize array ref
#                 $count_hr->{$chr}->{$start}->{'sum'} = 0;
#                 $count_hr->{$chr}->{$start}->{'nums'} = 0;
#             }
            $count_hr->{$chr}->{$start}->{'sum'} += $GC;
            $count_hr->{$chr}->{$start}->{'nums'}++;
        }
    }
    # now match GC in all and only the windows available in $hr
    $hr->{'gc'} = {};
    foreach $chr (keys(%{$hr->{'test'}})) {
        foreach my $start (keys(%{$count_hr->{$chr}})) {
            if (defined ($count_hr->{$chr}->{$start}) ) {
                $hr->{'gc'}->{$chr}->{$start} = 
                    $count_hr->{$chr}->{$start}->{'sum'} / 
                    $count_hr->{$chr}->{$start}->{'nums'};
            } else {
                $hr->{'gc'}->{$chr}->{$start} = 'NA';
            }
        }
    }
    return 1;
}

sub makeOutFileName {
    my ($testFile, $controlFile, $wSize) = @_;
    my($testFileName) = fileparse($testFile, qr/\.[^.]*/);
    my($controlFileName) = fileparse($controlFile, qr/\.[^.]*/);
    return $testFileName . "_vs_" . $controlFileName . "_" . $wSize . ".tab"; 
}

sub printHash {
    my ($hr) = @_;

    print "Chr\tPos\tTest\tNorm";
    print "\tGC" if ($hr->{'gc'}); 
    print "\n";
    my @chrs = sort chrSort keys(%{$hr->{'test'}});
    foreach my $chr (@chrs){
        my @s = sort{$a <=> $b} (keys(%{$hr->{'test'}->{$chr}}));
        foreach(@s) {
            print "$chr\t$_";
            print "\t$hr->{'test'}->{$chr}->{$_}";
            print "\t$hr->{'control'}->{$chr}->{$_}";
            if (defined ($hr->{'gc'})) {
                if (defined ($hr->{'gc'}->{$chr}->{$_})) {
                    print "\t$hr->{'gc'}->{$chr}->{$_}";
                } else {
                    print "\tNA";
                }
            }
            print "\n";

        }
    }
}

sub chrSort {
# sort words with a common alphabetic prefix follow by a number
# words with the same prefix will be numerically sorted according to the
# numeric part
    my $A = {}; 
    my $B = {}; 
    my $aa;
    my $bb;
    if ($a =~ /(\D*)(\d*)$/){
        $A->{'alfa'} = $1;
        $A->{'num'} = $2
    }
    if ($b =~ /(\D*)(\d*)$/){
        $B->{'alfa'} = $1;
        $B->{'num'} = $2;
    }
    if ($A->{'alfa'} eq $B->{'alfa'}){
        return $A->{'num'} <=> $B->{'num'};
    } elsif (defined($A->{'alfa'}) && defined($B->{'alfa'})) {
        return $A->{'alfa'} cmp $B->{'alfa'};
    } else {
        return $a cmp $b;
    }
}

sub sortNumFirst {
    if ($a =~ /\d+/ && $b =~ /\d+/){
        return $a <=> $b;
    } elsif ($a =~ /\d+/) {
        return -1;
    } elsif ($b =~ /\d+/) {
        return 1;
    } else {
        return $a cmp $b;
    }
}

sub equaliseAndFill {
    my ($hr1, $hr2, $wSize) = @_;   
    foreach my $chr (keys(%{$hr1})) {
        if (!defined($hr2->{$chr})) {
            die "$chr is missing from control sample!\n";
        }

    }
}

sub getHeaderFH {
    my ($file) = @_;
    my $fh;
    my $msg = "Impossible to open $file:";
    open ($fh, "samtools view -H '$file' |") || die "$msg $!\n";
    return $fh;
}

sub getFileHandle {
    my ($file) = @_;
    my $fh;
    my $msg = "Impossible to open $file:";
    if ($file =~ /\.gz$/) {
        open($fh, "gzip -dc '$file' |") || die "$msg $!\n";
    } elsif ($file =~ /\.bam$/) {
        # trying to open a bam file
        # first check samtools is in the path
        my $samtools = `which samtools`;
        chomp($samtools);
        if ($samtools){
            open($fh, "samtools view '$file' |") || die "$msg $!\n";
        } else {
            die "Trying to open a bam file, but `samtools` is not in \$PATH\n";
        }
    } else {
        open($fh, $file) || die "$msg $!\n";
    }
    
    return $fh;
}


sub setFilters {
    my (%par) = @_;
    my $filter_hr = {};
    my $blackList_r;
    my $threshold;
    if ($par{'blackListFile'}) {
        $blackList_r = loadBlackList($par{'blackListFile'});
    }
    if ($par{'qualityThreshold'}) {
        $threshold = $par{'qualityThreshold'};
    }
    $filter_hr->{'quality'} = $threshold;
    $filter_hr->{'blackList'} = $blackList_r;
    return $filter_hr;
}

sub isBlackListed {
    my ($seq_ar, $bl) = @_;
 
    if (!defined($bl->{$seq_ar->[2]})) {
        return 0;
    }
    for (my $n = 0; $n < @{$bl->{$seq_ar->[2]}->{'starts'}}; $n++) {
        if ($seq_ar->[3] >= $bl->{$seq_ar->[2]}->{'starts'}->[$n] &&
            $seq_ar->[3] < $bl->{$seq_ar->[2]}->{'ends'}->[$n]) {
            return 1;        
        }
    }   
    # never within a blacklist limit
    return 0;
}

sub thisSeqPassFilter {
    my ($seq_ar, $filters_hr) = @_;
    # as soon as it does not pass a filter, return zero.
    # Check faster filters first
    if (defined($filters_hr->{'quality'})) {
        if ($seq_ar->[4] < $filters_hr->{'quality'}) {
            return 0;
        }
    }
    if (defined($filters_hr->{'blackList'})) {
        if (isBlackListed($seq_ar, $filters_hr->{'blackList'})){
            return 0;
        }
    }
    # it was not filtered out, so it must be OK
    return 1;
}

sub line2pos {
    my ($line_ar, $filters_hr) = @_;
    if (thisSeqPassFilter($line_ar, $filters_hr)){
        return ($line_ar->[2], $line_ar->[3]);
    } else {
        return undef;
    } 

}

sub loadBlackList {
    my ($file) = @_; 
    my $hash_r = {}; 
    open (IN, $file) || die "Impossible to open $file: $!\n";
    my $firstLine = <IN>;
    while(<IN>) {
        chomp;
        my ($chr, $start, $end) = split(/\t/);
        if (!defined($hash_r->{$chr})) {
            $hash_r->{$chr}->{'starts'} = []; 
            $hash_r->{$chr}->{'ends'} = []; 
        }   
        push (@{$hash_r->{$chr}->{'starts'}}, $start);
        push (@{$hash_r->{$chr}->{'ends'}}, $end);
    }   
    return $hash_r;
}

sub _wcAdd {
    my ($what, $count_hr, $fh, $wSize) = @_;
    # this is a pseudoSub, uses variables from windowCount
    while (<$fh>){
        if (! /^\s*#/){
            chomp;
            my ($chr, $pos) = split(/\t/);
            if (defined ($chr)) {
                my $start = int($pos/$wSize) * $wSize + 1;
                $count_hr->{$what}->{$chr}->{$start}++;
#                 if (!defined($max_hr->{$chr}) || $start > $max_hr->{$chr}) {
#                     $max_hr->{$chr} = $start;
#                 } 
            }
        }
    }
}

sub windowCount {
    my ($wSize, $testFile, $controlFile, $filters_hr, $max_hr) = @_;
    my $count_hr = {};
    my $fh = getFileHandle($testFile);
    _wcAdd('test', $count_hr, $fh, $wSize);
    close $fh;
    $fh = getFileHandle($controlFile);
    _wcAdd('control', $count_hr, $fh, $wSize);
    close $fh;
    # Now check that all windows are the same in test and control
    foreach my $chr (keys(%{$count_hr->{'test'}})) {
        if (!defined($count_hr->{'control'}->{$chr})) {
            $count_hr->{'control'}->{$chr} = {};
        }

        for (my $n = 1; $n < $max_hr->{$chr}; $n = $n + $wSize) {
            if (!defined($count_hr->{'test'}->{$chr}->{$n})) {
                $count_hr->{'test'}->{$chr}->{$n} = 0;
            } 
            if (!defined($count_hr->{'control'}->{$chr}->{$n})) {
                $count_hr->{'control'}->{$chr}->{$n} = 0;
            } 
        }
    }
    return $count_hr;
}

sub checkSmaller {
    # recursive function. Checks if the previous window exists. If not creates
    # it and set it to zero
    my ($hr, $start, $wSize) = @_;
    my $prevStart = $start - $wSize;
    if ($prevStart < 1) {
        return 1;
    } elsif (defined($hr->{$prevStart})) {
        return 1;
    } else {
        $hr->{$prevStart} = 0;
        checkSmaller($hr, $prevStart, $wSize);
    }
}

sub setUsage {
    my %pars = @_;
    my $usage = qq|
$0 version $version

Copyright (C) 2010 Stefano Berri

Usage: perl $0 [Options] <testFile> <controlFile>

<testFile> Path to bam or (gzipped) sam file of test sample.
<controlFile> Path to bam or (gzipped) sam file of control sample.

Options:
    --window: size of window (in bp) to count reads. If this value is provided, 
        readNum will not be used. [$pars{'window'}]
    --gc_file. Path to a file with gc content as dowloaded from 
        UCSD. If a file is provided, GC content will be calculated  [$pars{'gc_file'}] 
    --readNum. Average number of reads in a window. Window size will
        be set accordingly. [$pars{'readNum'}]
    --genomeSize. Size of aploid genome. [$pars{'genomeSize'}]
    --qualityThreshold: sequences with MAPQ quality lower than qualityThreshold 
        will not be used (see bwa). Set to 0 for no filtering on quality score [$pars{'qualityThreshold'}]
    --tmpDir: path to temporary directory used for saving temporary files. [$pars{'tmpDir'}]
    --saveTest: a valid file name for a temporary file. If provided, the
        temporary file from the <testFile> will not be deleted and can be re-used
        with option --testTemp to save computing time [$pars{'saveTest'}]
    --saveControl: a valid file name for a temporary file. If provided, the
        temporary file from the <controlFile> will not be deleted and can
        be re-used with option --testTemp to save computing time. [$pars{'saveTest'}]
    --makeTempOnly: FLAG. If this flag is present, it only produce the temporary files
        and then exit. --saveTest and/or --saveControl must be present.
    --testTemp: FLAG. If this flag is present <testFile> is not a sam/bam
        but the *test* temporary file created using option --saveTest.
    --controlTemp:  FLAG. If this flag is present <controlFile> is not a
        sam/bam but the *control* temporary file created using option --saveControl.
\n|;
    return $usage;
# the following option is removed as not implemented yet.
# --blackListFile: Path to a tab delimited file with chromosome location
#             that should not considered. The file has a one line header and then 
#             tab delimited "chr\tstart\tend". [$pars{'blackListFile'}]
 
}

sub setDefaultPars {
    my %pars = (
        window => '',
        readNum => 30,
        gc_file => '',
        genomeSize => 3_000_000_000,
        qualityThreshold => 37,
        tmpDir => './',
        testTemp => undef,
        controlTemp => undef,
        saveTest => '',
        saveControl => '',
        makeTempOnly => undef,
        chrFile => '',
    );
    return %pars;
}


sub chrLength {
    my ($test, $control) = @_;
    my %chrLtest = getHead($test);
    my %chrLcontrol = getHead($control);
    if (hashCompare (\%chrLtest, \%chrLcontrol)){
        return (%chrLtest);
    } else {
        die "ERROR: Test and control refer to different reference genomes\n";
    }
}

sub hashCompare {
    my ($first_hr, $second_hr) = @_;
    if (lessOrEqualHash($first_hr, $second_hr) && lessOrEqualHash($second_hr, $first_hr)){
        return 1
    } else {
        return 0;
    }
}

sub lessOrEqualHash {
    my ($small_hr, $large_hr) = @_;
    foreach my $k (keys(%{$small_hr})){
        if (!defined($large_hr->{$k})){
            return 0;
        } elsif ($small_hr->{$k} ne $large_hr->{$k}){
            return 0;
        }
    }
    return 1;
}

sub getHead {
    my ($file) = @_;
    my %chrL;
    if (($file =~ /\.sam$/i) || ($file =~ /\.sam\.gz$/i)){
        my $fh = getFileHandle($file);
        while(<$fh>){
            chomp;
            my $thisLine = $_;
            if ($thisLine !~ /^\@/){
                return (%chrL);
            } else {
                addThisChrLength($thisLine, \%chrL);
            }
        }    
    } elsif ($file =~ /\.bam$/i) {
        my $fh = getHeaderFH($file);
        while(<$fh>){
            chomp;
            my $thisLine = $_;
            addThisChrLength($thisLine, \%chrL);
        }
        return (%chrL);
    } else {
        die "File $file is not valid\n";
    }
}

sub addThisChrLength {
    my ($thisLine, $chrL_hr) = @_;
    if ($thisLine =~ /\@SQ\tSN:(\S+)\tLN:(\d+)/){
        my $thisChr = $1;
        my $thisLength = $2;
        $chrL_hr->{$thisChr} = $2;
    } 
}

