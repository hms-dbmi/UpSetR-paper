input <- c(
  "MAQ"=144600,
  "FaSD"=16532, 
  "Bcftools"=283, 
  "GATK"=15160, 
  "MAQ&FaSD"=16323, 
  "MAQ&Bcftools"=636, 
  "Bcftools&GATK"=65435, 
  "FaSD&GATK"=33874, 
  "MAQ&FaSD&Bcftools"=114, 
  "MAQ&FaSD&GATK"=41858, 
  "MAQ&Bcftools&GATK"=4, 
  "FaSD&Bcftools&GATK"=6603, 
  "MAQ&FaSD&Bcftools&GATK"=8357
)

upset(fromExpression(input), order.by = "freq")
