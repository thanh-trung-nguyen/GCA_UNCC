# Comments:


# Example of use: 
# perl /home/thanhnguyen/UNCC_2012/programs/GCA/do_laplacetr.pl Sol_Z108.m laplacetr_z108_s ../Sim_UNCC_horn_hom_fre15/Sol_Z108.m 1200 0.003 51 51 5 8 7



$fileObject = $ARGV[0];
$fileOutput = $ARGV[1];
$fileHom = $ARGV[2];
$NoTimeSteps =   $ARGV[3];
$dt = $ARGV[4];
$Nx = $ARGV[5];
$Ny = $ARGV[6];
$s_min = $ARGV[7];
$s_max = $ARGV[8];
$Ns = $ARGV[9];



for ($i == 0; $i < $Ns; $i++)
{
	$s = $s_min + $i*($s_max - $s_min)/($Ns - 1);
	$Sname = sprintf('%1.2f', $s);
  
	`   /home/thanhnguyen/UNCC_2012/programs/GCA/laplace   $fileObject $NoTimeSteps $dt  $s $fileOutput$Sname.m $Nx $Ny $fileHom  `;


} 

