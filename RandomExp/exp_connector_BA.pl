#use strict;
#use warnings;

#print STDERR "Done.\n";
#exit(1);

@N = (1000,2000,5000,7500,10000,20000,30000,40000,50000); #頂点数
@NC = (50);  #初期クリークの頂点数
@AD = (3,5,10);  #頂点追加時に追加する辺の数
@Q = (25,50,100);  #アイテム数
@MQ = (5,10);   #頂点の持つ最大アイテム数
@SEED = (1,2,3,4,5); #シード値

@AD = (3,5);

$TLIM = 300;

foreach $numVer(@N){
   foreach $coreVer(@NC){
      foreach $addEdge(@AD){
         foreach $numItem(@Q){
            foreach $maxItem(@MQ){
               foreach $alg(("cooma","Tree", "Boley_et_al")){
			   $cmd = 'mkdir ./BA/result/'.$alg.'/'.$numVer.'_'.$coreVer.'_'.$addEdge.'_'.$numItem.'_'.$maxItem;
			   #print $cmd;
			   system($cmd);
			}
               foreach $seed(@SEED){
			my $cmd = ("python ./BA/gen_random_graph_BA.py mid.graph $seed $numVer $coreVer $addEdge \n");
			print $cmd;
			
			
			system($cmd);
			$cmd = "perl ./BA/Undirected_BA.pl mid.graph\n";
			print $cmd;
			system($cmd);

			$cmd = "python make_ptn.v03.py output.grh $numItem $maxItem output.ptn $seed\n";
			print $cmd;
			system($cmd);
			
			
			##$cmd = './copine/COPINE14 output.ptn output.grh 1 -t '.$TLIM.' > ./BA/result/copine/'.$numVer.'_'.$coreVer.'_'.$addEdge.'_'.$numItem.'_'.$maxItem.'/result_copine_'.$numVer.'_'.$coreVer.'_'.$addEdge.'_'.$numItem.'_'.$maxItem.'_'.$seed.'.txt';
			## print "$cmd\n";
			## system($cmd);

			$cmd = './BP/target/BoleyEtAl output.ptn output.grh 1 -tlim '.$TLIM.' > ./BA/result/Boley_et_al/'.$numVer.'_'.$coreVer.'_'.$addEdge.'_'.$numItem.'_'.$maxItem.'/result_Boley_'.$numVer.'_'.$coreVer.'_'.$addEdge.'_'.$numItem.'_'.$maxItem.'_'.$seed.'.txt';
			print "$cmd\n";
			system($cmd);

			$cmd = './FT/TREE2 output.ptn output.grh 1 -tlim '.$TLIM.' > ./BA/result/Tree/'.$numVer.'_'.$coreVer.'_'.$addEdge.'_'.$numItem.'_'.$maxItem.'/result_tree_'.$numVer.'_'.$coreVer.'_'.$addEdge.'_'.$numItem.'_'.$maxItem.'_'.$seed.'.txt';
			print "$cmd\n";
			system($cmd);

			$cmd = './cooma/A1B output.ptn output.grh 1 -tlim '.$TLIM.' > ./BA/result/cooma/'.$numVer.'_'.$coreVer.'_'.$addEdge.'_'.$numItem.'_'.$maxItem.'/result_cooma_'.$numVer.'_'.$coreVer.'_'.$addEdge.'_'.$numItem.'_'.$maxItem.'_'.$seed.'.txt';
			print "$cmd\n";
			system($cmd);

               }
            }
         }
      }
   }
}
