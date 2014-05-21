#!/hive/groups/gencode/local/bin/perl -w
use strict;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Registry;
use Getopt::Long;

# This parses the command line and sets the value of $all to 1 if the option
# was specified on the command line. 
my $all = 0;
GetOptions('all' => \$all);

my ($species, $chromFile, $outFile, @chroms);
$species = $ARGV[0];
$chromFile = $ARGV[1];
$outFile = $ARGV[2];

print "All is $all \n";

if ($#ARGV < 2) 
{
   print "Usage: getEnsembTxSeqsWithVersions.pl <species> <chroms list> <output file> \n";
   print "Chroms list can either be a list of chromosomes for which to retrieve transcript sequences or \'all\' to get all transcripts for an assembly.\n";
   print "Script either gets transcript sequences in FASTA format for specified chromosomes for the latest Ensembl build for a species or all Ensembl transcript sequences for that build.\n";
   exit 1;
}

# open chromosomes file and open output file for writing
if ($all != 1) 
{
   print "opening chroms file\n";
   open (CHROMS, $chromFile) or die "Cannot open $chromFile: $! \n";
   # Read in list of chromosomes
   @chroms = <CHROMS>;
   # Remove new lines and close file
   chomp(@chroms);
   close CHROMS;
}
open (OUT, ">$outFile") or die "Cannot create $outFile: $! \n";

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'useastdb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);
# get a slice adaptor for the human core database
my $slice_adaptor = $registry->get_adaptor( $species, 'Core', 'Slice' );
# get adaptor for transcripts
my $tr_adaptor = $registry->get_adaptor( $species, 'Core', 'Transcript');

# Print the transcript sequences in FASTA format
sub printTxSeqs
{
   my $tr_ref = shift(@_);
   while ( my $tr = shift@{$tr_ref} )
   {
      my $stableId = $tr->stable_id();
      my $version = $tr->version();
      my $seqObj = $tr->seq();
      my $txSeq = $seqObj->seq();
      print OUT ">$stableId.$version\n";
      print OUT "$txSeq\n";
   }
}

# Get all transcripts for each chromosome specified in a list
sub getAllTranscriptsByChrom
{
   my $transcripts;
   # Get a slice for each chromosome and process to get ids with versions. 
   for my $c (@chroms)
   {
      print "Chrom: $c\n";
      my $slice = $slice_adaptor->fetch_by_region('chromosome', $c);
      # Get all transcripts on this chromocome
      my $transcripts = $tr_adaptor->fetch_all_by_Slice($slice);
      # Print FASTA format sequences
      printTxSeqs($transcripts);
   }
}

# Gets all transcripts for an assembly
sub getAllTranscripts
{
   # Get reference to a list of the transcript objects
   my $transcripts = $tr_adaptor->fetch_all();
   printTxSeqs($transcripts);
}

# Get all the transcripts for specified organism if --all specified otherwise
# get transcripts for specified chromosomes and print FASTA sequences
if ($all == 1) 
{  
   print "Getting all transcript seqs \n";
   getAllTranscripts();
}
else 
{
   print "Getting transcript seqs for chroms \n";
   getAllTranscriptsByChrom();
}
