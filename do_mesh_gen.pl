# Comments:

# run the mesh_gen routine to generate the FEM mesh:



$N2 = $ARGV[0];


for ($i == 1; $i <= $N2; $i++)
{
	
	$inputfile = sprintf('%s%d%s', '../parameter_vis_obj',$i,'.dat');
	print $inputfile; print "\n";
  
	`   /home/tnguy152/UNCC_2012/programs/GCA/mesh_gen $inputfile `;


} 

