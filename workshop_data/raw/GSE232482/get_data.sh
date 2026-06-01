wget --recursive --no-parent -nd ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232482/suppl/

tar -xvf GSE232482_RAW.tar

ls | grep "Human" | xargs rm

for smp in `ls *_singlecell.tar.gz`; do
  smp_name=`basename ${smp} _filtered_peak_bc_matrix_fragments_singlecell.tar.gz`
  mkdir -p ${smp_name}
  mv ${smp} ${smp_name}
  cd ${smp_name}
  tar -xzvf ${smp}
  cd /igm/home/cnh008/ref/GSE232482
done

for smp in `ls *_matrix.tar.gz`; do
  smp_name=`basename ${smp} _filtered_feature_bc_matrix.tar.gz`
  mkdir -p ${smp_name}
  mv ${smp} ${smp_name}
  cd ${smp_name}
  tar --strip-components=10 -xzvf ${smp}
  cd /igm/home/cnh008/ref/GSE232482
done

# god bless america there are so many dang nested folders. modified to add the --strip-components to remove the !! 10 !! subdirectories