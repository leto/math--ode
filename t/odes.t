use Test::More tests => 3;
use File::Spec;
use lib File::Spec->catfile("..","lib");
use Math::ODE;
use Data::Dumper;
my ($o, $eps, @vals, $res);

# analytic solution is y(x) = 2 e^{-x}
$o = new Math::ODE ( file => 'data',
			step => 0.1,
			initial => [2], 
			DE => [ \&DE1 ], 
			t0 => 0,
			tf => 5 );
$o->evolve;
$eps = $o->error;	
@vals =  $o->values_at(3);
$res = abs( $vals[0]  - 2*exp(-3) );
ok( $res < $eps, "Constant Coefficient 1st order solved correctly, res=$res, eps=$eps"); 


# analytic solution is y(x) = -1/(z-1)
$o = new Math::ODE ( file => 'data',
			step => 0.01, 
			initial => [1], 
			DE => [ \&DE2 ], 
			t0 => 0,
			tf => 0.5 );
$o->evolve;
$eps = $o->error;	
@vals =  $o->values_at(0.4);
$res = abs( $vals[0]  + 1/(0.4-1) );
ok( $res < $eps, "Nonlinear 1st order solved correctly, res=$res, eps=$eps"); 


# analytic solution is y(x) = exp(-x^2)
$o = new Math::ODE ( file => 'data',
			step => 0.1, 
			initial => [1], 
			DE => [ \&DE3 ], 
			t0 => 0,
			tf => 5 );
$o->evolve;
$eps = $o->error;	
@vals =  $o->values_at(3);
$res = abs( $vals[0]  + 2*3*exp(-3**2) );
ok( $res < $eps, " 1st order nonhomogeneous solved correctly, res=$res, eps=$eps"); 

sub DE1 { my ($t,$y) = @_; return -$y->[0]; }
sub DE2 { my ($t,$y) = @_; return $y->[0] ** 2; }
sub DE3 { my ($t,$y) = @_; return - 2 * $y->[0] * exp ( - $y->[0] ** 2 ); }
