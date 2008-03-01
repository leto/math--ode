use Test::More tests => 25;
use File::Spec;
use lib File::Spec->catfile("..","lib");
use Math::ODE;
use Data::Dumper;
use strict;
use warnings;

my ($o, $ok, $sol, $error, @vals, $res);
my @steps = qw(0.1 0.09 0.08 0.05 0.01 0.005);

######################################
# y'   = - y 
# y(0) = 2
######################################
# analytic solution is y(x) = 2 e^{-x}
######################################
for my $step (@steps) {
    $o = new Math::ODE ( 
            step => $step,
            initial => [2], 
            ODE => [ \&DE1 ], 
            t0 => 0,
            tf => 5 );
    if( $o->evolve ){
        $sol = sub { my $t=shift; 2 * exp(-$t) };
        $error = $o->max_error([$sol]);
        ok( $error < $o->error , "y'=y, y(0)=2, y(x) = 2 e^{-x}\nstep=$step, max error=$error, expected error=" .$o->error); 
    } else {
        ok( 0 );
    }
}

######################################
# y' = y^2
# y(0) = 1 
######################################
# analytic solution is y(x) = -1/(x-1)
######################################
for my $step (@steps) {
    $o = new Math::ODE ( 
			step => $step, 
			initial => [1], 
			ODE => [ \&DE2 ], 
			t0 => 0,
			tf => 0.5 );
    $o->evolve;

    $sol = sub { my $t=shift; -1/($t-1) };
    $error = $o->max_error([$sol]);

    ok( $error < $o->error ,  "y'=y^2, y(0)=1, y(x)=-1/(x-1)\n step=$step, max error=$error, expected error=" . $o->error); 
}

######################################
# y' = - 2 x exp(-x^2)
# y(0) = 1
######################################
# analytic solution is y(x) = exp(-x^2)
######################################
sub DEgauss { my ($t,$y) = @_; - 2 * $t * exp ( - $t ** 2 ); }
for my $step (@steps) {
    $o = new Math::ODE ( 
			step => $step,
			initial => [1], 
			ODE => [ \&DEgauss ], 
			t0 => 0,
			tf => 5 );
    $o->evolve;

    $sol = sub { my $t=shift;  exp(-$t**2) };

    $error = $o->max_error([$sol]);

    ok( $error < $o->error , "y'=-2 x exp(-x^2), y(0) = 1, y(x) = exp(-x^2)\n step=$step, max error=$error, expected error=" . $o->error); 
}
##############
##############

sub DE1 { my ($t,$y) = @_; -$y->[0]; }
sub DE2 { my ($t,$y) = @_; $y->[0] ** 2; }

SKIP: {

skip "stuff", 1;

sub DE3 { my ($t,$y) = @_; $y->[0] ** 3; }
sub DE4 { my ($t,$y) = @_; $y->[0] ** 4; }
sub DE5 { my ($t,$y) = @_; $y->[0] ** 5; }
sub DE6 { my ($t,$y) = @_; $y->[0] ** 6; }

for my $step (@steps) {
        for my $n ( 3 .. 4 ){
            my $subname = "DE$n";
            $o = new Math::ODE (    
                        step => $step, 
                        initial => [1], 
                        ODE => [\&$subname],
                        t0 => 0,
                        tf => 0.5 );
            $o->evolve;
            my $sol = sub { my $x = shift; (1-$n)*($x+1/(1-$n)) ** (1/($n-1)) };
            $error = $o->max_error([$sol]);
            ok( $error < $o->error , "$n-th order nonlinear: step=$step, max error=$error, expected error=" . $o->error); 

        }
}

}
# y' = x^-1, y(1) = 1 
#  Solution: y = ln(x)

sub DElog { my ($t,$y) = @_; return 1/$t; }
for my $step (@steps) {
    $o = new Math::ODE (
            step => $step,
            initial => [0], 
            ODE => [ \&DElog ], 
            t0 => 1,
            tf => 10 );
    if ( $o->evolve ) {
        $sol = sub { log(shift) };
        $error = $o->max_error([$sol]);
	    ok( $error < $o->error , "y'=1/x, y(1)=1, y(x)=log(x)\n step=$step, max error=$error, expected error=" .$o->error); 
    } else {
	    ok( 0, 'numerical badness');
    }

}

