use Template;
 
# some useful options (see below for full list)
my $config = {
    minia => '/search/path/minia3',  # or list ref
    abyss2 => '/search/path/abyss2',  # or list ref
    disco => '/search/path/DiscoVarDenovo',  # or list ref
    FM => '/search/path/fastmin-sg',  # or list ref
    IM => '/search/path/intervalmiss',  # or list ref
    LIGER => '/search/path/liger',  # or list ref
    AUXPS => '/search/path/scripts',  # or list ref
    seqtk=>'/search/path/seqtk',
};
 
# create Template object
my $template = Template->new($config);
 
# define template variables for replacement
my $vars = {
    minia => undef,  # or list ref
    abyss2 => undef,  # or list ref
    disco => undef,  # or list ref
    FM => '/Users/adigenova/Git/fastmin-sg/fastmin-sg',  # or list ref
    IM => '/Users/adigenova/Git/intervalmis/intervalmis',  # or list ref
    LIGER => '/Users/adigenova/Git/liger/liger',  # or list ref
    AUXPS => undef,  # or list ref
    seqtk=>'/Users/adigenova/Git/seqtk/seqtk',
};
 
# specify input filename, or file handle, text reference, etc.
my $input = 'Global.tt';
 
# process input template, substituting variables
$template->process($input, $vars)
    || die $template->error();
