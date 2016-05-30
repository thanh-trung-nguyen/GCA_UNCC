# Comments:

# run the GCA for the full data set:



for ($i == 0; $i <= 27; $i++)
{
	
	$object = sprintf('%s%d', 'object_',$i);
	chdir $object;
	print $object, "\n";
  
	` ../../GCA/gca_totalwave parameter_forprob.dat  ../parameter_inversion_ub_nonmetal.dat `;
	
	chdir;

} 

