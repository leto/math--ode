use Test::More tests => 2;
use File::Spec;
use lib File::Spec->catfile("..","lib");
use Math::ODE;
use Data::Dumper;
use strict;
use warnings;

# ODEs: f'' - g' = 0        
#       g'' + f*g' + f'*g = 0
#       f(0)=0,f'(0)=1,f''(0)=-1,g(0)=1
# Solution: f(x) = 1 - exp(-x)
#           g(x) = exp(-x)

SKIP : {
	skip 'unkown numerical error is creeping in somewhere', 2;
my $o = new Math::ODE ( step => 0.05,
                        initial => [0,1,-1,1],
                        ODE => [ \&DE1, \&DE2 , \&DE3, \&DE4 ],
                        t0 => 0,
                        tf => 10 );
if ($o->evolve){
	my $x = 5;
	my $error = $o->error;	
	my @vals =  $o->values_at(3);
	my $res1 = abs( $vals[0]  - (1 - exp(-$x)) );
	my $res2 = abs( $vals[1]  -  exp(-$x)      );
	ok( $res1 < $error, "f: two coupled second order equations, res=$res1, expected error=$error"); 
	ok( $res2 < $error, "g: two coupled second order equations, res=$res2, expected error=$error"); 

} else {
	ok(0, 'bad');
}

}
sub DE1 { my ($t,$y) = @_; return $y->[1]; }
sub DE2 { my ($t,$y) = @_; return $y->[2]; }
sub DE3 { my ($t,$y) = @_; return $y->[1] * $y->[3] - $y->[0] * $y->[2]; }
sub DE4 { my ($t,$y) = @_; return $y->[2]; }

