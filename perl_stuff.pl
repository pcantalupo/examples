
# glob
my @files = glob ("*.bam");
print join ("\n", @files), "\n";

# command line
# $0 is script name
# @ARGV contains command line arguments
print join (" ", $0, @ARGV),"\n";
