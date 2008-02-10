use Test::More tests => 2;
use File::Spec;
use lib File::Spec->catfile("..","lib");
use Math::ODE;
use Data::Dumper;

# analytic solution is y(x) = 5 x
my $o = new Math::ODE ( file => 'data',
			step => 0.1,
			initial => [0], 
			DE => [ \&DE1 ], 
			t0 => 0,
			tf => 1 );
ok( ref $o eq 'Math::ODE', 'new returns proper object' );
$o->evolve;
my $eps = $o->{step} ** 4;	# because Math::ODE implements a 4th order Runge-Kutta method

my $s = sprintf("%0.12f", 0.5);
my @vals =  $o->values_at( $s );
my $res = abs($vals[0] - 2.5);
ok( $res < $eps, "Constant Coefficient Equation solved correctly, res=$res"); 


sub DE1 { my ($t,$y) = @_; return 5; }

