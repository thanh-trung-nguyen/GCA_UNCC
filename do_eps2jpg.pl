# Comments:

# convert EPS figures to JPG figures;



for ($i == 1; $i <= 18; $i++)
{
	
	$inputfile = sprintf('%s%d%s','obj',$i,'_3D.eps');
	$outputfile = sprintf('%s%d%s','obj',$i,'_3D.jpg');
  
	`   convert $inputfile $outputfile`;

	$inputfile = sprintf('%s%d%s','obj',$i,'_xy.eps');
	$outputfile = sprintf('%s%d%s','obj',$i,'_xy.jpg');
  
	`   convert $inputfile $outputfile`;

	$inputfile = sprintf('%s%d%s','obj',$i,'_yz.eps');
	$outputfile = sprintf('%s%d%s','obj',$i,'_yz.jpg');
  
	`   convert $inputfile $outputfile`;


} 

