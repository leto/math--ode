use File::Spec;
use lib File::Spec->catfile("..","lib");
use Test::More tests => 1;
use Math::ODE;

my $o = new Math::ODE ( file => 'data',
			step => 0.1,
			initial => [0], 
			DE => [ \&DE1 ], 
			t0 => 0,
			tf => 1 );
ok( ref $o eq 'Math::ODE', 'new returns proper object' );
$o->evolve;
print join ", ", @{ $o->{values}{0} }; 
print join ", ", @{ $o->{values}{0.1} }; 
print join ", ", @{ $o->{values}{0.5} }; 
print join ", ", @{ $o->{values}{1.0} }; 

my $stuff = $o->file;
print "stuff=$stuff\n";

sub DE1 { my ($t,$y) = @_; return 5; }

