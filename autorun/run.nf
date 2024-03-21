
process create_GLASS_OD_metadata__step_001 {

  input:
    // file "../glass-od-clinical-database/glass-od-clinical-database.db"
    // file "scripts/load_GLASS-OD_metadata__step_001.R"
  
  output:
    path "cache/glass_od.metadata.patients__step_001.Rds"
    path "cache/glass_od.metadata.resections__step_001.Rds"
    path "cache/glass_od.metadata.array_samples__step_001.Rds"

  
  """
  export wd=`pwd`;
  echo \$wd;
  
  cd $projectDir;
  Rscript --vanilla scripts/load_GLASS-OD_metadata__step_001.R
  
  
  mkdir -p \$wd/cache ;
  cp $projectDir/cache/glass_od.metadata.patients__step_001.Rds      \$wd/cache/glass_od.metadata.patients__step_001.Rds ;
  cp $projectDir/cache/glass_od.metadata.resections__step_001.Rds    \$wd/cache/glass_od.metadata.resections__step_001.Rds ;
  cp $projectDir/cache/glass_od.metadata.array_samples__step_001.Rds \$wd/cache/glass_od.metadata.array_samples__step_001.Rds ;
  """
}

workflow {
    create_GLASS_OD_metadata__step_001()
}
