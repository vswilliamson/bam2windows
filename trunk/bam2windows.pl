# Prepare input file for CNAnorm discrete normalization from a pair of sam/bam
# files. If you use bam format, samtools must be in your $PATH or you need to
# pass its path with --samtools-path
# Type:
# perl <nameOfThisScript>
# to see usage
#
#
# Copyright (C) 2010 Stefano Berri <sberri@illumina.com> or <s.berri@gmail.com>
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
# created 25/08/2010
# last update 09/01/2014
my $version = "0.3.9";



use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;

my %par = setDefaultPars();
my $usage = setUsage(%par);
my $verboseVersion = setVerboseVersion($version);

GetOptions (\%par, 'version|v', 'verbose|V', 'window|w=i', 'gc_file|gc=s',
    'samtools-path=s', 'readNum|r=i', 'genomeSize=i', 'qualityThreshold|q=i',
    'tmpDir|d=s', 'testTemp|tt', 'controlTemp|ct', 'testSorted|ts',
    'controlSorted|cs', 'makeTempOnly|t', 'saveTest|st=s', 'saveControl|sc=s',
    'chrFile=s', 'help|h');

if ($par{'version'}){print $verboseVersion; exit 0}
if ($par{'help'}){print $usage; exit 0}

%par = parameterCheck(%par);

if (@ARGV == 0 || @ARGV > 2){print $usage; exit 1}
if (@ARGV == 1 && !defined($par{'makeTempOnly'})){print $usage; exit 2}

# create filter object
if ($par{'verbose'}){
    warn "initialiasing filters...\n"
}
my $filters_hr = setFilters(%par);

my $testFile = shift;
my ($countTestFile, $readNumTest) = sam2simple($testFile, $par{'tmpDir'}, 
    $par{'testTemp'}, $filters_hr, $par{'saveTest'}, $par{'verbose'});
# exit if (!defined($controlFile) && $par{'makeTempOnly'});

my $controlFile = shift;
if (!defined($controlFile) && $par{'makeTempOnly'}){
    # we only wanted to create a temp file from one sample. Fine. Exit
    exit 0;
}

my ($countControlFile, $readNumControl) = sam2simple($controlFile, $par{'tmpDir'},
    $par{'controlTemp'}, $filters_hr, $par{'saveControl'}, $par{'verbose'});

# if we don't need only the temp files...
if (!defined($par{'makeTempOnly'})){
    
    # load default chromosome length
    # my %chrLength = chrLenght($par{'chrFile'});
    my %chrLength = chrLength($testFile, $controlFile, $par{'testTemp'}, $par{'controlTemp'});
    addGenomesize(\%par, \%chrLength);

# how many reads in the smallest of the two files?
    my $NumReads = min($readNumTest, $readNumControl);
    if (!$par{'window'}) {
        # set window size accordingly
        $par{'window'} = getWindowSize($NumReads, $par{'genomeSize'},
            $par{'readNum'});
    }


    # produce counts
    if($par{'verbose'}) {
        warn "Counting reads over window $par{'window'} bp wide...\n";
    }
    my $hr = windowCount($par{'window'}, $countTestFile, $countControlFile, 
        $par{'testSorted'}, $par{'controlSorted'}, \%chrLength, $par{'verbose'});
    
    if ($par{'gc_file'} ne '') {
        if ($par{'verbose'}){
            warn "calculating GC content...\n";
        }
        addGC_content($hr, $par{'gc_file'}, $par{'window'});        
    }
    printHash($hr, $par{'verbose'}); 

}


# rename or delete temporary file
if (! renameOrDelete($countTestFile, $par{'saveTest'}, 
    $par{'testTemp'}, $par{'verbose'})) {
    warn "Impossible to rename or delete test temporary file\n";
}

if ($controlFile){
    if (! renameOrDelete($countControlFile, $par{'saveControl'}, 
        $par{'controlTemp'}, $par{'verbose'})) {
        warn "Impossible to rename or delete control temporary file\n";
    }
}

if ($par{'verbose'}){
    warn "Done!\n";
}



#### SUBS ####

sub addGenomesize {
    my ($par_hr, $chrL_hr) = @_;
    my $totSize = 0;
   
    foreach (keys %{$chrL_hr}){
        $totSize += $chrL_hr->{$_};
    }
    $par_hr->{'genomeSize'} = $totSize;
}

sub renameOrDelete {
    my ($file, $saveFile, $tempFile, $verbose) = @_;
    my $success;
    if ($tempFile){
        $success = 1;
    } elsif ($saveFile ne ''){
        $success = rename $file, $saveFile;
    } else {
        if ($verbose){
            warn "removing file $file...\n";
        }
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
    my ($file, $dir, $tmpF, $filters_hr, $saveFile, $verbose) = @_;
    my $numOfGoodReads = 0;
    
    if ($tmpF){
        open (TMP, $file) || die "Impossible to open temporary file $file\n";
        if($verbose){
            warn "$file is a previously produced temporary file, recovering total numer of reads...\n";
        }
        my $firstLine = <TMP>;
        if ($firstLine =~ /^# 0*(\d+)/){
            seek (TMP, $1, 0);
        } else {
            die "file $tmpF is not a valid temporary file"
        }
        
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
        
        if($verbose){
            warn "$file is a sam/bam file. Converting, filtering and recovering total numer of reads...\n";
        }
        my ($fileName) = fileparse($file);
        my $fh = getFileHandle($file, $par{'samtools-path'});
        $dir =~ s/\/$//;
        my $tmpFile;
        if ($saveFile){
            $tmpFile = $dir . "/" . $saveFile;
        } else {
            $tmpFile = $dir . "/" . createTmpFileName($fileName);
        }
        # warn "creating temporary file $tmpFile...\n";
        open (OUT, ">$tmpFile") || die "Impossible to create temporary file $tmpFile: $!\n";
        # prepare space for pointer. Will be substitute at the end
        print OUT "# 000000000000\n";
        while (<$fh>){
            if (/^\@.*/){
                print OUT "#$_"; 
            } else {
                my @F = split(/\t/);
                my ($chr, $pos) = line2pos(\@F, $filters_hr);
                if (defined ($chr)) {
                    print OUT "$chr\t$pos\n";
                    $numOfGoodReads++;
                }
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
    my ($hr, $verbose) = @_;
    if ($verbose){
        warn "printing file...\n";
    }
    print "Chr\tPos\tTest\tNorm";
    print "\tGC" if ($hr->{'gc'}); 
    print "\n";
    my @chrs = sort chrSort keys(%{$hr->{'test'}});
    foreach my $chr (@chrs){
        my @s = sort{$a <=> $b} (keys(%{$hr->{'test'}->{$chr}}));
        foreach(@s) {
            print "$chr\t$_";
            if (defined($hr->{'test'}->{$chr}->{$_})){
                print "\t$hr->{'test'}->{$chr}->{$_}";
            } else {
                print "\t0";
            }
            if (defined($hr->{'control'}->{$chr}->{$_})){
                print "\t$hr->{'control'}->{$chr}->{$_}";
            } else {
                print "\t0";
            }
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
    my ($file, $samtools) = @_;
    my $fh;
    my $msg = "Impossible to open $file:";
    open ($fh, "$samtools view -H '$file' |") || die "$msg $!\n";
    return $fh;
}

sub getFileHandle {
    my ($file, $samtools) = @_;
    my $fh;
    my $msg = "Impossible to open $file:";
    if ($file =~ /\.gz$/) {
        open($fh, "gzip -dc '$file' |") || die "$msg $!\n";
    } elsif ($file =~ /\.bam$/) {
        # trying to open a bam file
        # first check samtools is in the path
        open($fh, "$samtools view -h '$file' |") || die "$msg $!\n";
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

sub wcAdd {
    my ($what, $count_hr, $fh, $wSize, $sorted) = @_;
    my $start = 1;
    my $end = $wSize;
    my $thisWinCount = 0;
    my $prevChr = '';
    while (<$fh>){
        if (! /^#/) {
            chomp;
            my ($chr, $pos) = split(/\t/);
            if ($sorted){
                if (($chr eq $prevChr && $pos <= $end) || $prevChr eq ''){ 
                    $thisWinCount++;
                } else {
                    $count_hr->{$what}->{$prevChr}->{$start} = $thisWinCount;
                    $start = int($pos/$wSize) * $wSize + 1;
                    $end = $start + $wSize - 1;
                    $thisWinCount = 1;
                }
                $prevChr = $chr;
            } else { 
                $start = int(($pos - 1)/$wSize) * $wSize + 1;
                $count_hr->{$what}->{$chr}->{$start}++;
            }
        }
    }
    # flush last field
    if ($sorted){
        $count_hr->{$what}->{$prevChr}->{$start} = $thisWinCount;
    }
}

sub makeEmptyHash {
    my ($wSize, $max_hr, $what_ar) = @_;
    my $hr = {};
    foreach my $chr (keys(%{$max_hr})) {
        for (my $n = 1; $n <= $max_hr->{$chr}; $n = $n + $wSize){
            foreach my $w (@{$what_ar}){
                $hr->{$w}->{$chr}->{$n} = 0;
            }
        }
    }
    return ($hr);
}

sub windowCount {
    my ($wSize, $testFile, $controlFile, $tSorted, $cSorted, $max_hr, $verbose) = @_;
    
    my $count_hr = makeEmptyHash($wSize, $max_hr, ['test', 'control']);
    my $fh = getFileHandle($testFile, $par{'samtools-path'});
    if ($verbose){
        warn "  counting reads per window for 'test' file $testFile...\n";
    }
    wcAdd('test', $count_hr, $fh, $wSize, $tSorted);
    close $fh;
    $fh = getFileHandle($controlFile, $par{'samtools-path'});
    if ($verbose){
        warn "  counting reads per window for 'control' file $controlFile...\n";
    }
    wcAdd('control', $count_hr, $fh, $wSize, $cSorted);
    close $fh;
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


sub chrLength {
    my ($test, $control, $testTmp, $controlTmp) = @_;
    my %chrLtest = getHead($test, $testTmp);
    my %chrLcontrol = getHead($control, $controlTmp);
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
    my ($file, $isTmp) = @_;
    my %chrL;
    if ($isTmp){
        my $fh = getFileHandle($file, $par{'samtools-path'});
        my $null = <$fh>; # this should always exist in temp file
        while(defined(my $line = <$fh>) && $null){
            chomp($line);
            if ($line =~ /^#(\@SQ.*)/){
                addThisChrLength($1, \%chrL);
            } else {
                $null = 0; # to exit the while loop  
            }
        }
        return (%chrL);
    } elsif (($file =~ /\.sam$/i) || ($file =~ /\.sam\.gz$/i)){
        my $fh = getFileHandle($file, $par{'samtools-path'});
        while(<$fh>){
            chomp;
            my $thisLine = $_;
            if ($thisLine !~ /^\@SQ/){
                return (%chrL);
            } else {
                addThisChrLength($thisLine, \%chrL);
            }
        }    
    } elsif ($file =~ /\.bam$/i) {
        my $fh = getHeaderFH($file, $par{'samtools-path'});
        while(<$fh>){
            chomp;
            my $thisLine = $_;
            addThisChrLength($thisLine, \%chrL);
        }
        return (%chrL);
    } else {
        die "File $file is not a valid sam/bam file\n";
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
    if ($par{'genomeSize'}){
        warn "Option 'genomeSize' is deprecated, the genome size will be calculated from\n" .
            "sum of chromosome length as describe in sam/bam header\n";
    }
    # check about samtools
    if ($par{'samtools-path'}){
        if (! -f $par{'samtools-path'}){
            die "The specified path to samtools ($par{'samtools-path'}) is not a valid file\n";
        }
        if (! -x $par{'samtools-path'}){
            die "The specified path to samtools ($par{'samtools-path'}) is not executable\n";
        } 
    } else { # path was not defined, should be in $PATH
        my $samtools = `which samtools`;
        chomp($samtools);
        if ($samtools){
            $par{'samtools-path'} = $samtools;
        } else {
            die "samtools is not in \$PATH and no alternative was specified with --samtools-path\n";
        }
    }
    return %par;
}

sub setVerboseVersion {
    my ($version) = @_;
    my $verboseV = qq|
bam2windows version $version.
Copyright (C) 2010 Stefano Berri

This is free software.  You may redistribute copies of it under the terms of
the GNU General Public License <http://www.gnu.org/licenses/gpl.html>.
There is NO WARRANTY, to the extent permitted by law.
|;
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
    -v, --version: FLAG. Print version and exit (status 0)
    -V, --verbose: FLAG. Warn about progress
    -h, --help: show this help and exit (status 0)
    -w, --window: size of window (in bp) to count reads. If this value is provided, 
        readNum will not be used. [$pars{'window'}]
    -gc, --gc_file. Path to a file with gc content as dowloaded from 
        UCSD. If a file is provided, GC content will be calculated  [$pars{'gc_file'}] 
    --samtools-path: Path to samtools. If not defined, it has to be in \$PATH or
        the program will abort [$pars{'samtools-path'}]
    -r, --readNum. Average number of reads in a window. Window size will
        be set accordingly. [$pars{'readNum'}]
    -q, --qualityThreshold: sequences with MAPQ quality lower than qualityThreshold 
        will not be used (see bwa). Set to 0 for no filtering on quality score [$pars{'qualityThreshold'}]
    -d, --tmpDir: path to temporary directory used for saving temporary files. [$pars{'tmpDir'}]
    -st, --saveTest: a valid file name for a temporary file. If provided, the
        temporary file from the <testFile> will not be deleted and can be re-used
        with option --testTemp to save computing time [$pars{'saveTest'}]
    -sc, --saveControl: a valid file name for a temporary file. If provided, the
        temporary file from the <controlFile> will not be deleted and can
        be re-used with option --testTemp to save computing time. [$pars{'saveTest'}]
    -t, --makeTempOnly: FLAG. If this flag is present, it only produce the temporary files
        and then exit. --saveTest and/or --saveControl must be present.
    -tt, --testTemp: FLAG. If this flag is present <testFile> is not a sam/bam
        but the *test* temporary file created using option --saveTest.
    -ct, --controlTemp:  FLAG. If this flag is present <controlFile> is not a
        sam/bam but the *control* temporary file created using option --saveControl.
    -ts, --testSorted: FLAG. If the test input file is sorted, set this flag to save
        computational time
    -cs, --controlSorted: FLAG. If the control input file is sorted, set this flag to save
        computational time
    
    DEPRECATED options
    --genomeSize. Size of aploid genome. [$pars{'genomeSize'}]
    --chrFile. This option is now ignored
\n|;
    return $usage;
}

sub setDefaultPars {
    my %pars = (
        window => '',
        verbose => 0,
        help => 0,
        readNum => 30,
        gc_file => '',
        'samtools-path' => '', 
        genomeSize => '',
        qualityThreshold => 37,
        tmpDir => './',
        testTemp => undef,
        controlTemp => undef,
        testSorted => undef,
        controlSorted => undef,
        saveTest => '',
        saveControl => '',
        makeTempOnly => undef,
        chrFile => '',
    );
    return %pars;
}


