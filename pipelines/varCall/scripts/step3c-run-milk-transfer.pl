#!/usr/bin/perl -w

use strict;
use Time::HiRes qw(usleep nanosleep);
use File::Path;
use File::Basename;
use FindBin;
use lib $FindBin::Bin;
use hyunlib qw(%hszchrs initRef forkExecWait);
use gcconfig;

&initRef($ref);

my @chrs = @ARGV;

if ( $#chrs < 0 ) {
    die __FILE__ . ": Usage: [command] [list of chromosomes separated by space]\n";
}


my $modelchr = "1";
if($chrs[0] eq "X"){$modelchr = "X";}
my $modelchrsz = $hszchrs{$modelchr}->[3];

my $modelf = "$out/svm/$modelchr\_1_$modelchrsz\_milk_svm.svm.model";
my $milkDir = "$out/aux/milk";
my $svmDir = "$out/svm";
my $pasteDir = "$out/paste";

my @ignores = qw(AC AN AF GC GN HWEAF HWDGF);
my @includes = ();
my $checkNA = 1;
my $cutoff = -0.5;

## hard filter criteria 
## HWE p-value p<1e-6
## 

unless ( -e "$svmDir" ) {
    mkdir("$svmDir") || die __FILE__ . ": Cannot create directory $svmDir\n";
}

my @names = ();
my $ncols = 0;
my %hIgnores = map { $_ => 1 } @ignores;
my %hIncludes = ();
my @includeKeys = ();
my @includeDefaultValues = ();
for(my $i=0; $i < @includes; ++$i) {
	my ( $key, $defaultVal ) = split(/,/,$includes[$i]);
	$hIncludes{$key} = $i;
	push(@names,$key);
	push(@includeKeys,$key);
	push(@includeDefaultValues,defined($defaultVal) ? $defaultVal : 0);
}
my $nIncludes = $#includes + 1;

my @discords = (0,1,1,0,0,1,1,0,1,0,0,1,0,0,0,1,0,0,1,0,1,1,0,0,1,1,0);
my @dupdisc = (0,1,1,1,0,1,1,1,0);	

foreach my $chr (@chrs) {
    my $szchr = $hszchrs{$chr}->[3];
    my $prefix = "$svmDir/$chr\_1\_$szchr\_milk_trans";;
    my $base_prefix = "$chr\_1\_$szchr\_milk_transfer";
    # simplify output file name for the sake of snakemake
    my $simplified_prefix = "$chr\_milk_transfer";

    print STDOUT __FILE__ . ": Writing the feature information for libsvm..\n";
    open(RAW,">$prefix.raw") || die "Cannot open $prefix.raw for writing\n";
    open(SITE,">$prefix.site") || die "Cannot open $prefix.raw for writing\n";
    open(IN,"zcat $pasteDir/$chr\_1\_$szchr\_paste.sites.vcf.gz|") || die "Cannot open file\n";

    my @hdrs = ();
    while(<IN>) {
	if ( /^#/ ) {
	    push(@hdrs,$_);
	    next;
	}
	my ($chrom,$pos,$id,$ref,$alt,$qual,$filt,$info) = split(/[ \t\r\n]+/);
	$info .= ";QUAL=$qual" unless ( $qual eq "." );
	my @infos = split(/;/,$info);
	my @values = ();
	my $k = 0;
	for(my $j=0; $j < @includeKeys; ++$j) {
	    push(@values,$includeDefaultValues[$j]);
	    ++$k;
	}	    
	
	for(my $j=0; $j < @infos; ++$j) {
	    my ($key,$val) = split(/=/,$infos[$j]);
	    next if ( defined($hIgnores{$key}) );
	    next unless defined($val); ## skip keys without any values
	    if ( defined($hIncludes{$key}) ) {
			if ( !($checkNA) || ($val =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) ) {
			    $values[$hIncludes{$key}] = $val; # set value if given
		}
	    }
	    else {
		if ( !($checkNA) || ( $val =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) ) {
		    push(@values,$val);
		}
		else {
		    push(@values,0);
		}
		
		if ( $ncols == 0 ) {
		    push(@names,$key);
		}
		else {
		    die "Cannot recognize $key in $infos[$j], supposed to be $names[$j] at $j, $chrom:$pos:$ref:$alt $info\n" unless ($names[$k] eq $key );
		}
		++$k;
	    }
	}
	if ( $ncols == 0 ) {
	    $ncols = $#names+1;
	    print STDOUT __FILE__ . ": Recording following $ncols features : @names\nThe info field was $info\n@infos\n";
	    die if ( $ncols == 0 );
	}
	elsif ( $ncols != $#values+1 ) {
	    die "Number of columns are not identical at $chrom:$pos:$ref:$alt\n";
	}
	print SITE join("\t",$chrom,$pos,$id,$ref,$alt,$qual,$filt,$info)."\n";	    
	print RAW join("\t",@values)."\n";
    }
    close IN;
    close RAW;
    close SITE;

    print STDOUT __FILE__ . ": Performing quantile normalization of features..\n";    
    my $cmd = "$invNorm --in $prefix.raw --out $prefix.norm";
    &forkExecWait($cmd);

    open(SITE,"$prefix.site") || die "Cannot open $prefix.sites\n";
    open(NORM,"$prefix.norm") || die "Cannot open $prefix.norm\n";
    open(FTR,">$prefix.feature") || die "Cannot open $prefix.feature\n";

    my ($npos,$nneg,$noth) = (0,0,0);
    while(<SITE>) {
	my ($chrom,$pos,$id,$ref,$alt) = split;
	my @z = split(/[ \t\r\n]+/,<NORM>);
	my $ln = "";
	for(my $i=0; $i < @names; ++$i) {
	    $ln .= " ".($i+1).":$z[$i]";
	}
	print FTR "0 $ln\n";
	++$noth;
    }
    close FTR;
    close SITE;

    print STDOUT __FILE__ . ": Applying trained SVM model on $noth variants\n";

    $cmd = "$svmclassify $prefix.feature $modelf $prefix.svm.pred";
    &forkExecWait($cmd);

    print STDOUT __FILE__ . ": Writing filtered site VCF files with SVM scores..\n";
    open(SITE,"$prefix.site") || die "Cannot open $prefix.sites\n";    
    open(PRED,"$prefix.svm.pred") || die "Cannot open $prefix.svm.pred file\n";
    open(OUT,"| $bgzip -c > $svmDir/${simplified_prefix}.sites.vcf.gz") || die "Cannot open $svmDir/$base_prefix.sites.vcf.gz\n";
    open(OUTU,"| $bgzip -c > $prefix.uniq.sites.vcf.gz") || die "Cannot open $prefix.uniq.sites.vcf.gz\n";

    my $milkVcf = "$milkDir/$chr\_1\_$szchr\_milk.sites.vcf.gz";
    unless ( -e $milkVcf ) { die "File " . $milkVcf . " not found!"  }
    open(MILK,"zcat $milkVcf | cut -f 1-8| grep -v ^#|") || die "Cannot open file\n";    

    splice(@hdrs,$#hdrs,0,"##INFO=<ID=SVM,Number=1,Type=Float,Description=\"Milk-SVM score for variant quality, passing -0.5 or greater\">\n");
    splice(@hdrs,$#hdrs,0,"##FILTER=<ID=SVM,Description=\"Variant failed SVM filter\">\n");
    splice(@hdrs,$#hdrs,0,"##INFO=<ID=OVERLAP,Number=.,Type=String,Description=\"Overlap with possible variant types\">\n");
    splice(@hdrs,$#hdrs,0,"##FILTER=<ID=DISC,Description=\"Mendelian or duplicate genotype discordance is high (3 or more)\"\n");
    splice(@hdrs,$#hdrs,0,"##FILTER=<ID=EXHET,Description=\"Excess heterozygosity with HWE p-value < 1e-6\"\n");     
    print OUT join("",@hdrs);
    splice(@hdrs,$#hdrs,0,"##FILTER=<ID=CONFLICT,Description=\"Variant position conflicts with other variants with higher SVM score\"\n");
    print OUTU join("",@hdrs);    
    
    my %hbest = ();
    
    while(<SITE>) {
	my @F = split;
	my $pred = <PRED>;

	chomp($pred);

	my $milk = <MILK>;
	chomp($milk);

	my @M = split(/[\t\r\n ]+/,$milk);
	while( defined($M[4]) && ( $M[4] =~ /,/ ) ) {
	    $milk = <MILK>;
	    chomp($milk);
	    @M = split(/[\t\r\n ]+/,$milk);	    
	}

	my @tcnts = split(/,/,$1) if ( $M[7] =~ /TRIO_CONC_THRES=([^;]+)/ );
	my @dcnts = split(/,/,$1) if ( $M[7] =~ /DUP_CONC_THRES=([^;]+)/ );	
	my ($ndisc,$nconc,$ddisc,$dconc) = (0,0,0,0);
	for(my $i=1; $i < 27; ++$i) {
	    if ( $tcnts[$i] > 0 ) {
		if ( $discords[$i] == 0 ) { $nconc += $tcnts[$i] }
		else { $ndisc += $tcnts[$i]; }
	    }
	}
	for(my $i=1; $i < 9; ++$i) {
	    if ( $dcnts[$i] > 0 ) {
		if ( $dupdisc[$i] == 0 ) { $dconc += $dcnts[$i] }
		else { $ddisc += $dcnts[$i]; }		
	    }
	}

	my $hweslp = $1 if ( $F[7] =~ /HWE_SLP=([^;]+);/ );

	my @filts = ();
	if ( ( $pred < $cutoff ) || ( $F[7] =~ /AC=0;/ ) ) {
	    push(@filts,"SVM");
	}
	if ( ( $ddisc > 2 ) || ( $ndisc > 2 ) ) {
	    push(@filts,"DISC");
	}
	if ( $hweslp > 13.81551 ) {
	    push(@filts,"EXHET");	    
	}
	    
	$F[6] = ($#filts < 0) ? "PASS" : join(';',@filts);
	$F[7] .= ";SVM=$pred";
	$F[7] .= ";OVERLAP=".join(",",split(/;/,$M[6])) if ( ( $M[6] ne "." ) && ( $M[6] ne "PASS" ) );
	$F[7] =~ s/QUAL=([^;]+);//;
	
	print OUT join("\t",@F)."\n";

	#next if ( length($F[3]) + length($F[4]) > 2 );
	if ( defined($hbest{$F[1]}) ) {
	    if ( $hbest{$F[1]}->[2] < $pred ) {
		$hbest{$F[1]} = [$F[3],$F[4],$pred];		
	    }
	}
	else {
	    $hbest{$F[1]} = [$F[3],$F[4],$pred];
	}
    }
    close OUT;
    close PRED;
    close SITE;

    my $ns = 19879;
    my $acbreaks = join(",",1,2,3,int($ns*0.002),int($ns*0.02),int($ns*0.2));
# Show the parameters for vcfSummary
   print "vcfsummary  $vcfsummary\n";
   print "ref  $ref\n";
   print "dpsnp  $dbsnp\n";
   print "hapmapvcf  $hapmapvcf\n";
   print "bgzip   $bgzip\n";
   print "tabix  $tabix\n";
   print "chr  $chr\n";
   print "out $prefix\n";
 
    print STDOUT __FILE__ . ": Producing VCF summary for SNPs..\n";

    open(IN,"zcat $svmDir/$simplified_prefix.sites.vcf.gz|") || die "Cannot open file\n";
    open(OUT1, "| perl $vcfsummary --ref $ref --db $dbsnp --FNRvcf $hapmapvcf --bgzip $bgzip --tabix $tabix --chr $chr --info AC --info-breaks $acbreaks > $prefix.sites.summary_ac") || die "Cannot open file\n";
    open(OUT2, "| perl $vcfsummary --ref $ref --db $dbsnp --FNRvcf $hapmapvcf --bgzip $bgzip --tabix $tabix --chr $chr > $prefix.sites.summary") || die "Cannot open file\n";
    open(OUT3, "| perl $vcfsummary2 --ref $ref --db $dbsnp --FNRvcf $hapmapvcf --bgzip $bgzip --tabix $tabix --chr $chr > $prefix.sites.summary_v2") || die "Cannot open file\n";    
    while(<IN>) {
	next if ( /^#/ );
	my @F = split;
	print OUT3 $_;
	if ( ( defined($hbest{$F[1]}) ) && ( $hbest{$F[1]}->[0] eq $F[3] ) && ( $hbest{$F[1]}->[1] eq $F[4] ) ) {
	    print OUTU $_;
	}
	else {
	    print OUTU join("\t",@F[0..5],$F[6].";CONFLICT",$F[7])."\n";
	    next;
	}
	#next if ( length($F[3]) + length($F[4]) > 2 );	
	#print OUT1 $_;
	print OUT2 $_;
    }
    close IN;
    close OUT1;
    close OUT2;
    close OUT3;
    close OUTU;
}

