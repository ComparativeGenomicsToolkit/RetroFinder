#!/hive/groups/gencode/local/bin/perl -w
use strict;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Registry;

my ($chromFile, @chroms);
$chromFile = $ARGV[0];
open (CHROMS, $chromFile) or die "Cannot open $chromFile: $! \n";
# Read in list of chromosomes
@chroms = <CHROMS>;
# Remove new lines and close file
chomp(@chroms);
close CHROMS;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

# get a slice adaptor for the human core database
my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
# get adaptor for transcripts
my $tr_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Transcript');

# Get a slice for each chromosome and process to get ids with versions. 
for my $c (@chroms)
{
   print "Chrom: $c\n";
   my $slice = $slice_adaptor->fetch_by_region('chromosome', $c);
   # Get all transcripts on this chromocome
   my $transcripts = $tr_adaptor->fetch_all_by_Slice($slice);
   while ( my $tr = shift @{$transcripts} ) 
   {
      my $dbId = $tr->dbID();
      my $stableId = $tr->stable_id();
      my $version = $tr->version();
      my $seqObj = $tr->seq();
      my $txSeq = $seqObj->seq();
      print "Chrom: $c Transcript $stableId [$dbId] Version: $version \n";
      print "Chrom: $c Transcript $stableId.$version \n";
      #print ">$stableId.$version\n";
      #print "$txSeq\n";
   }
}

# copy of list in reference
# my @transcriptsList = @{$transcripts_ref};

