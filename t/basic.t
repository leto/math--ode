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
my $vals =  $o->values(0.5); 

ok( abs( $vals->[0] - 2.5) < $eps, "Constant Coefficient Equation solved correctly, eps=$eps"); 


sub DE1 { my ($t,$y) = @_; return 5; }

