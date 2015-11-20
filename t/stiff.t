use strict;
use warnings;
use File::Spec::Functions;
use Test::More;
use File::Spec;
use lib File::Spec->catfile("..","lib");
use Math::ODE;
use Data::Dumper;
my $file = catfile(qw/t data/);

# This is an extremely stiff ODE can our normal integrator barfs on
# More details in examples/stiff_de1

my $A=499999.5;
my $B=500000.5;
my $de1 = sub { my ($t,$y) = @_; -$A*$y->[0] + $B*$y->[1] };
my $de2 = sub { my ($t,$y) = @_; -$B*$y->[0] - $A*$y->[1] };

my $o = new Math::ODE(
        file    => 'data',
        step    => 0.5,
        initial => [0,2],
        ODE     => [ $de1, $de2],
        t0      => 0,
        tf      => 10 );
$o->evolve;

done_testing;

END { unlink $file };
