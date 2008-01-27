use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Math::ODE',
    license             => 'perl',
    dist_author         => 'Jonathan Leto <jonathan@leto.net>',
    dist_version_from   => 'lib/Math/ODE.pm',
    dist_abstract	=> 'Solve N-th Order Differential Equations (Initial Value Problems)',
    build_requires => {
        'Test::More' => 0,
    },
    add_to_cleanup      => [ 'Math::ODE-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
