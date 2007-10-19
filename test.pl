# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..4\n"; }
END {print "not ok 1\n" unless $loaded;}
use Math::ODE;
$loaded = 1;
print "ok 1\n";

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):


my $o = new Math::ODE ( file => 'data',
			verbose => 1,
			step => 0.1,
			initial => [0,1],     #, -1, 1 ],
			DE => [ \&DE1, \&DE2 ], #, \&DE3, \&DE4 ],
			t0 => 0,
			tf => 10 );
print "ok 2\n";
$o->evolve();
print "ok 3\n";
$o->file("data");
print "ok 4\n";


#########
=head
sub DE1 { my ($t,$y) = @_; return $y->[1]; }
sub DE2 { my ($t,$y) = @_; return $y->[2]; }
sub DE3 { my ($t,$y) = @_; return $y->[1] * $y->[3] - $y->[0] * $y->[2]; }
sub DE4 { my ($t,$y) = @_; return $y->[2]; }
=cut

# initial [0,1]
# plot 'data' using 1:2, sin(x)
#=head  y'' + y = 0
sub DE1 { my ($t,$y) = @_; return $y->[1]; }
sub DE2 { my ($t,$y) = @_; return -$y->[0]; }
#=cut

# initial [0,1]
# plot 'data' using 1:2, x*exp(-x)
=head y'' + 2y' + y = 0
sub DE2 { my ($t,$y) = @_; return -2 * $y->[1] - $y->[0]; }
sub DE1 { my ($t,$y) = @_; return $y->[1]; }
=cut
