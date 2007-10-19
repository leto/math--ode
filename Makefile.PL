use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'		=> 'Math::ODE',
    'VERSION'		=> '0.03', 
    ($] >= 5.005 ?    
      (ABSTRACT => 'Solve N-th Order Ordinary Differential Equations',
       AUTHOR     => 'Jonathan Leto <jonathan@leto.net>') : ()),
    'PREREQ_PM'		=> { 'Data::Dumper' => 0 }, 
);
