use Test::More tests => 2;
use File::Spec;
use lib File::Spec->catfile("..","lib");
use Math::ODE;
use Data::Dumper;

# analytic solution is y(x) = 2 e^{-x}
my $o = new Math::ODE ( file => 'data',
			step => 0.1,
			initial => [2], 
			DE => [ \&DE1 ], 
			t0 => 0,
			tf => 5 );
$o->evolve;
print Dumper [ $o ];
my $eps = $o->{step} ** 4;	# because Math::ODE implements a 4th order Runge-Kutta method
my @vals =  @{ $o->{values}{3} }; 
my $res = abs( $vals[0]  - 1 );
ok( $res < $eps, "Constant Coefficient Equation solved correctly, res=$res"); 


sub DE1 { my ($t,$y) = @_; return -$y->[0]; }

