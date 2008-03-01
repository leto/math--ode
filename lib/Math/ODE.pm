package Math::ODE;
use strict;

require 5.005;
use Data::Dumper;
use Carp;
my $VERSION = '0.04';

$Data::Dumper::Varname = "y";
$Data::Dumper::Indent = 0;

sub evolve {
	my $self = shift;
	my ($F,$h,$t,$y,$file) = map{ $self->{$_} } qw(ODE step t0 initial file);
	my $delim = $self->{csv} ? ',' : ($self->{delim} || $self->{delimeter} || " ");
	
	if( defined $file ){
        	open(FD, ">$file") or croak "$file: $!";
	}

    while ( $t < $self->{tf} ){
        # use Runge Kutta 4th order to step from $t to $t + $h
        $y = _RK4($self,$t,$y);

        warn "Exiting RK4 with t=$t ," . Dumper($y) . "\n" if( $self->{verbose} > 1 );

        for my $i ( 0 .. $self->{N}-1 ){
            # check for under/over flow
            next unless $y->[$i] =~ qr/nan|infinity/i;
            warn "Bailing out, over/under flow at t=$t,y->[$i] = $y->[$i]" if $self->{verbose};
            return undef;
        }
        $t += $h;

        if( defined $file ){
            my $str = join $delim,  map { sprintf "%0.12f", $_ } ($t, @$y);
		    chop $str;
            print FD "$str\n";
	    }
    }
	close FD if defined $file;
    return 42;
}

sub _RK4 {
        # $t = dependent variable
        # $y = $N - vector of independent variables
        # $h = step size
        # $F = arrayref of coderefs of the equations to solve
        my $self = shift;
        my ($t, $y) = @_;
        my $F = $self->{ODE};
        my $h = $self->{step};

        ## w vectors hold constants for equations
        ## each $q holds a modified $y vector to feed to the next
        ## for loop ( ie $y + $w1/2 , etc ... )
        my (@w1,@w2,@w3,@w4,$q,$i);

        for $i ( 0 .. $self->{N}-1 ){ $w1[$i]  = $h * &{ $F->[$i] }($t,$y);           }
        for $i ( 0 .. $self->{N}-1 ){ $q->[$i] = $y->[$i] + 0.5*$w1[$i];              }

        for $i ( 0 .. $self->{N}-1 ){ $w2[$i]  = $h * &{ $F->[$i] }($t + 0.5*$h,$q);  }
        for $i ( 0 .. $self->{N}-1 ){ $q->[$i] = $y->[$i] + 0.5*$w2[$i];              }

        for $i ( 0 .. $self->{N}-1 ){ $w3[$i]  = $h * &{ $F->[$i] }($t + 0.5*$h,$q);  }
        for $i ( 0 .. $self->{N}-1 ){ $q->[$i] = $y->[$i] + $w3[$i];                  }

        for $i ( 0 .. $self->{N}-1 ){ $w4[$i]  = $h * &{ $F->[$i] }($t + $h,$q);      }


        for $i ( 0 .. $self->{N}-1 ){ $y->[$i] += ( $w1[$i] + 2 * $w2[$i] + 2 * $w3[$i] + $w4[$i])/6; }

	    $self->_store_values( $t + $h, $y );
	
        return $y;
}
sub _store_values {
	my ($self,$t, $y) = @_;
	return unless  $self->{keep_values};
	my $s = sprintf '%0.12f', $t ; 
	push @{ $self->{values}{$s} }, @$y;
}
sub values_at {
	my ($self,$t, %args) = @_;
    if ($self->{keep_values}){
	    return @{ $self->{values}{sprintf('%0.12f',$t)} };
    } else {
        warn "Values were not kept because keep_values was set to 0";
        return;
    }
}
# because Math::ODE implements a 4th order Runge-Kutta method
sub error {  $_[0]->{step} ** 4 }

sub _init {
	my ($self,%args) = @_;

	# defaults
	$self->{keep_values} = 1;
	$self->{verbose}     = 1;
	$self->{step}        = 0.1; 
	$self->{csv}         = 0;
	$self->{N}           = scalar( @{ $args{ODE} } ) || 1;

	@$self{keys %args} = values %args;
	$self->{values} = {};

    if( $self->{N} != scalar(  @{ $args{initial } }) ){
                croak "Must have same number of initial conditions as equations!";
    }
	if( $self->{step} <= 0  ){
		croak "Stepsize must be positive!";
	}
	if( $self->{t0} >= $self->{tf} ){
		croak "\$self->t0 must be less than \$self->tf!";
	}
    return $self;
}
sub max_error {
    my ($self, $sols, $eps) = @_;
    $eps ||= 1e-8;

    my $maxres = 0;
    my (@vals, $res);

    for my $pt ( sort keys %{$self->{values}} ){
        my $k = 0;
        for my $sol ( @$sols ) {
            @vals =  $self->values_at($pt);
            #next unless defined $vals[0] && &$sol($pt);
            $res = abs( $vals[$k]  - &$sol($pt) );
            $maxres = $res if ($res > $maxres);
            #print "pt=$pt, res=$res\n" if ($res > $self->error && debug() );
            $k++;
        }
    }
    $maxres;
}
sub debug { 0 }

sub new {
	my $class = shift;
	my $self = {};
	bless($self, $class);
	$self->_init(@_);
}
42;
__END__

=head1 NAME

Math::ODE - Solve N-th Order Ordinary Differential Equations 

=head1 DESCRIPTION

This module allows you to solve N-th Order Ordinary Differential Equations with
as little pain as possible.  Currently, only IVP's (initial value problems) are
supported, but native support for BVP's (boundary value problems) may be added
in the future. To solve N-th order equations, you must first turn it into a 
system of N first order equations, as in MATLAB.

=head1 SYNOPSIS

	use Math::ODE;
	# create new object that stores data in a file 
	# and solve the given equation(s) on the interval
	# [0,10], with initial condition y(t0) = 0
	my $o = new Math::ODE ( file => '/home/user/data',
                        step => 0.1,
                        initial => [0], 
                        ODE => [ \&DE1 ], 
                        t0 => 1,
                        tf => 10 );
	$o->evolve();
	# solve the equation y' = 1/$t
	# $t is the independent variable, a scalar
	# $y is the dependent variable, an array reference
	sub DE1 { my ($t,$y) = @_; return 1/$t; }

=over 2

=item *

C<$o-E<gt>evolve()>

Evolves the equations from C<$o-E<gt>t0> to C<$o-E<gt>tf> using a 4th Order Classic Runge-Kutta method.

	# Example 1: Solve  y'' + y = 0, y(0) = 0, y'(0) = 1
	#            Solution: y = sin(x)

	use Math::ODE;
	my $o = new Math::ODE ( file => 'data',
				verbose => 1,
				step => 0.1,
				initial => [0,1], 
				ODE => [ \&DE0, \&DE1 ], 
				t0 => 0,
				tf => 10 );
	$o->evolve;

	# to plot data in gnuplot: plot 'data' using 1:2, sin(x)

	# y'' + y = 0
	sub DE0 { my ($t,$y) = @_; return $y->[1]; }
	sub DE1 { my ($t,$y) = @_; return -$y->[0]; }

To turn y'' + y = 0 into a system, we will imagine a vector with two
components, called y. Now let the first component y0 = y and the second
component y1 = y0' .  Now rewrite the equation in terms of these variables. It
will be y1' + y0 = 0. Solving for y1' we get y1' = -y0. Now the vector y' has
compenents y0' = y1 and y1' = -y0 . These are the equations we put in our DE*
sub's, except we use an arrayref. The DE0 sub corresponds to our y0' component
and the DE1 sub corresponds to our y1' component. This can be *ahem* easily 
generalized to any N-th Oder equation or system of equations. An example of
a system of second order equations is given in example/de4. Please look in
there for other examples as well.

The C<initial> arrayref corresponds to the initial conditions of each of the
components of the dependent variable vector. C<[0,1]> means that the first
component at C<$o-E<gt>t0> (which happens to be 0) will have the value 0 and
the second component at C<$o-E<gt>t0> (which is 10) will have the value 1.

=item *

C<$o-E<gt>ODE($F)>

Sets the equations to be solved to C<$F>. C<$F> must be an arrayref
of coderefs and the same length as C<$F>, or Bad Things May Occur (tm). 

=item *

C<$o-E<gt>initial($arrayref)>

Set the initial conditions to C<$arrayref>. Returns an arrayref
of initial conditions if no arguments are given. Must be the same
size as C<$o-E<gt>ODE> or autovivisection will ensue.

=item *

C<$o-E<gt>step($s)>

This will set the step size (length of the segments of the
interval) to C<$s>, which must be positive. The default is 0.01.
Returns the step size if no arguments are given.

=item *

C<$o-E<gt>t0($a)>

This will set the initial (leftmost) point on the interval to C<$a>.
Returns the initial point on the interval if no arguments are given.

=item *

C<$o-E<gt>tf($b)>

This will set the last (rightmost) point on the interval to C<$b>.
Returns the end point on the interval if no arguments are given.

=item *

C<$o-E<gt>file($somefile)>

Save data in $somefile. Returns the file in which the data is being saved if no
arguments are given. If no file is specified, data is printed to STDOUT.  The
data file will have $N+1 columns (where $N is the number of equations to
solve), and can be fed directly to gnuplot. The first column is the independent
variable, and the remaining are the first through nth components of the
dependent vector. Examples of graphing the data file are in the example/
directory of the source distribution.

=item *

C<$o-E<gt>verbose($number)>

Sets the verbosity of debugging output. The default of 1 will currently only
cause C<evolve()> to C<warn> you when it gets C<NaN> or C<Inf>, instead of 
silently returning C<undef>, when verbosity is set to 0. Setting the verbosity
to 0 is useful if you are trying to write code to solve boundary value problems
with the shooting method. You will be guessing initial conditions, and don't feel
like getting a warning on every guess that blows up (which will be many.) Setting
the verbosity to 2 will cause a message like the following:

	Exiting RK4 with t=9.9 ,$y1 = ['-0.544013766248772833','-0.839075464413064726'];

to be printed on every increment of the independent variable C<$t>. These are the values
that the 4th Order Runge-Kutta returned for the current value of C<$t>.



=back

=head1 AUTHOR

Jonathan Leto <jonathan@leto.net>

=head1 SEE ALSO

Boyce, DiPrima "Elementary Differential Equations" 5th Ed.

Orwant, Hietaniemi, Macdonald "Mastering Algorithms with Perl" Ch. 16.

=head1 COPYRIGHT

Copyright (c) 2001-2008 by Jonathan Leto.  All rights reserved.

=head1 LICENSE AGREEMENT

This package is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=cut
