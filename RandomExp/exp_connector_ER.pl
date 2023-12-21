#use strict;
#use warnings;


@N = (1000,2000,5000,7500,10000,20000,30000,40000,50000); #頂点数
@D = (3,5,10,15);           #辺の生成比率（論文のd）
@Q = (25,50,100);         #アイテム数: 10 is cancelled.
@MQ = (5,10);         #頂点の持つ最大アイテム数
@SEED = (1,2,3,4,5);        #シード値


$TLIM = 300;

foreach $numVer(@N){
   foreach $ratioEdge(@D){
         foreach $numItem(@Q){
	     foreach $maxItem(@MQ){
		 if($maxItem > $numItem){
		     print STDERR "warning: q=$numItem, dq=$maxItem is skipped.\n";
		     next;
		 }
               foreach $alg(("cooma","Tree", "Boley_et_al")){
			   $cmd = 'mkdir ./ER/result/'.$alg.'/'.$numVer.'_'.$ratioEdge.'_'.$numItem.'_'.$maxItem;
			   #print $cmd;
			   system($cmd);
		}
               my $ratio = $ratioEdge / ($numVer - 1); 
               foreach $seed(@SEED){
			my $cmd = ("python ./ER/gen_random_graph_ER_20230221.py mid.graph $seed $numVer $ratio \n");
			print $cmd;
			
			
			system($cmd);
			$cmd = "perl ./ER/Undirected_ER.pl mid.graph\n";
			print $cmd;
			system($cmd);
			$cmd = "python make_ptn.v03.py output.grh $numItem $maxItem output.ptn $seed\n";
			print $cmd;
			system($cmd);
			
			



			##$cmd = './copine/COPINE14 output.ptn output.grh 1 -t '.$TLIM.' > ./ER/result/copine/'.$numVer.'_'.$ratioEdge.'_'.$numItem.'_'.$maxItem.'/result_copine_'.$numVer.'_'.$ratioEdge.'_'.$numItem.'_'.$maxItem.'_'.$seed.'.txt';
			##print "$cmd\n";
			##system($cmd);

			$cmd = './BP/target/BoleyEtAl output.ptn output.grh 1 -tlim '.$TLIM.' > ./ER/result/Boley_et_al/'.$numVer.'_'.$ratioEdge.'_'.$numItem.'_'.$maxItem.'/result_Boley_'.$numVer.'_'.$ratioEdge.'_'.$numItem.'_'.$maxItem.'_'.$seed.'.txt';
			print "$cmd\n";
			system($cmd);

			$cmd = './FT/TREE2 output.ptn output.grh 1 -tlim '.$TLIM.' > ./ER/result/Tree/'.$numVer.'_'.$ratioEdge.'_'.$numItem.'_'.$maxItem.'/result_tree_'.$numVer.'_'.$ratioEdge.'_'.$numItem.'_'.$maxItem.'_'.$seed.'.txt';
			print "$cmd\n";
			system($cmd);

			$cmd = './cooma/A1B output.ptn output.grh 1 -tlim '.$TLIM.' > ./ER/result/cooma/'.$numVer.'_'.$ratioEdge.'_'.$numItem.'_'.$maxItem.'/result_cooma_'.$numVer.'_'.$ratioEdge.'_'.$numItem.'_'.$maxItem.'_'.$seed.'.txt';
			print "$cmd\n";
			system($cmd);

            
            }
         }
      }
   }
}
