library(fitdistrplus)

# Input: read counts file
# 1. fit the allele frequency by a beta distribution
# 2. fit the depth of coverage by a negative binomial
# 3. fit the loci expression level by a gamma distribution
# Return: the distribution parameters for the three fits
getNBparam <- function(file){
  RC = read.table(file)
  colnames(RC) = c("genome_RC_A", "genome_RC_B",
                   "transcript_RC_A","transcript_RC_B")
  
  # 1. Genome coverage is fitted by a negative binomial distribution
  fit_Gcoverage = fitdist(RC$genome_RC_A + RC$genome_RC_B, 
                          method = "mme", distr="nbinom")
  plot(fit_Gcoverage)
  GenomeCov_NB_theta = fit_Gcoverage$estimate[1]
  GenomeCov_NB_mu = fit_Gcoverage$estimate[2]
  
  # 2. Genome variant frequency distribution is fitted by a beta distribution
  fit_freq = fitdist( RC$genome_RC_A / (RC$genome_RC_A + RC$genome_RC_B), 
                      method = "mme", distr="beta")
  VarFreq_Beta_shape1 = fit_freq$estimate[1]
  VarFreq_Beta_shape2 = fit_freq$estimate[2]
  plot(fit_freq)
  
  # 3. Loci expression level is fitted by a gamma distribituion
  fit_Tcoverage = fitdist(RC$transcript_RC_A + RC$transcript_RC_B, 
                          distr="gamma", method = "mme")
  Expression_Gamma_shape = fit_Tcoverage$estimate[1]
  Expression_Gamma_rate = fit_Tcoverage$estimate[2]
  plot(fit_Tcoverage)
  
  
  return (as.numeric(as.character(c(GenomeCov_NB_theta, GenomeCov_NB_mu, 
                                    VarFreq_Beta_shape1, VarFreq_Beta_shape2, 
                                    Expression_Gamma_shape, Expression_Gamma_rate ))))
}

# One read count file per sample
files = list.files ("RC")
files = paste("RC/",files,sep="")
parameters = matrix(unlist(lapply(files,getNBparam)), ncol=6, byrow = TRUE)
colnames(parameters) = c("GenomeCov_NB_theta", "GenomeCov_NB_mu", 
                         "VarFreq_Beta_shape1", "VarFreq_Beta_shape2", 
                         "Expression_Gamma_shape", "Expression_Gamma_rate")
rownames(parameters) = c("S155","S158", "S178","S206","S208","S209","S210")
write.table(x = parameters, "genomic_distributions_parameters.txt", sep = "\t", quote = FALSE)
print (parameters)

library(ggplot2)
p = ggplot(data = as.data.frame(parameters), 
           aes( x = GenomeCov_NB_mu, y = GenomeCov_NB_theta))+
  geom_point()+
  geom_smooth(method = lm)+
  xlab(expression(mu))+
  ylab(expression(theta))+
  theme_bw()
print(p)

# Predict the NB theta  with mu by a linear model
fitParam = summary(lm(parameters[,1]~parameters[,2]))
print (fitParam)

mu2theta <- function (mu){
  theta = fitParam$coefficients[1] + fitParam$coefficients[2] * mu
  return (theta)
}

test_RC_Noise <- function (genome_coverage, theta, fbeta1, fbeta2, gamshape, gamrate, transcript_sd){
  # simulate random allele frequency given the beta distrib parameters (shape1 and shape2)
  variant_frequency = rbeta(n = 1, shape1 = fbeta1, shape2 = fbeta2)
  genomic_expect_cov_A = round( variant_frequency * genome_coverage , 0)
  genomic_expect_cov_B = genome_coverage - genomic_expect_cov_A
  
  ## simulate genomic read counts for allele A and B 
  # given expected coverages and the NB distrib parameters
  # first find theta
  thetaA = round(mu2theta(genomic_expect_cov_A))
  thetaB = round(mu2theta(genomic_expect_cov_B))
  # then simulate
  genomic_RC_A = rnbinom (n = 1, size = thetaA, mu = genomic_expect_cov_A) 
  genomic_RC_B = rnbinom (n = 1, size = thetaB, mu = genomic_expect_cov_B)
  # replace NA by 0
  if (is.na(genomic_RC_A)) genomic_RC_A = 0 
  if (is.na(genomic_RC_B)) genomic_RC_B = 0
  
  # simulate allele expression level given the gamma distrib parameters (shape and rate)
  gene_expression_Level = round (rgamma(1, gamshape, gamrate), 0)
  transcriptomic_expect_Cov_A = round(variant_frequency * gene_expression_Level , 0)
  transcriptomic_expect_Cov_B = gene_expression_Level - transcriptomic_expect_Cov_A
  
  # simulate transcriptomic read counts for allele A and B
  # given a poisson distribution
  transcript_RC_A =  rpois(n = 1, transcriptomic_expect_Cov_A )
  transcript_RC_B =  rpois(n = 1, transcriptomic_expect_Cov_B )
  if (transcript_RC_A<0) transcript_RC_A = 0 
  if (transcript_RC_B<0) transcript_RC_B = 0
  
  # Apply Fisher exact test 
  # only if the gene is expressed and the loci is heterozygous
  observed_expression = transcript_RC_A + transcript_RC_B
  if( (observed_expression>0 ) & genomic_RC_A !=0 & genomic_RC_B !=0 ){
    Fisher_pvalue = fisher.test( matrix(c(genomic_RC_A, 
                                          genomic_RC_B, 
                                          transcript_RC_A, 
                                          transcript_RC_B),
                                        nrow=2))$p.value
  }
  else{
    Fisher_pvalue = "NA"
  }
  return(as.numeric(as.character(c(Fisher_pvalue, genomic_RC_A , genomic_RC_B, 
                                   transcript_RC_A, transcript_RC_B))))
}

n_test = 50000
results = c()
# To compare to psADE found in Oithona %ADE for q_value <0.05
oithona = c(0.21, 0.14, 0.33, 1.15, 0.27, 0.03, 0.43)
names(oithona) = rownames(parameters)

for (sample in rownames(parameters)){
  FisherTest = matrix(c(rep(0, 5*n_test)), ncol=5)
  
  for (j in 1:n_test){
    FisherTest[j,] = test_RC_Noise(round(parameters[sample,"GenomeCov_NB_mu"],0), thetaMu,
                                   parameters[sample,"VarFreq_Beta_shape1"], 
                                   parameters[sample,"VarFreq_Beta_shape2"],
                                   parameters[sample,"Expression_Gamma_shape"],
                                   parameters[sample,"Expression_Gamma_rate"],
                                   parameters[sample,"Expression_variation_sd"]
    )
  }
  # Writesimulation output
  dir.create("simulated_data")
  colnames(FisherTest) = c("Fisher_pvalue", "genomic_RC_A" , "genomic_RC_B", 
                           "transcript_RC_A", "transcript_RC_B")
  write.table(FisherTest, paste("simulated_data/",sample,".txt", sep = ""))
  
  # Check if the genomic read count simulation is OK
  print (sample)
  genome_simulated_coverage = FisherTest[,2] + FisherTest[,3]
  # expected distribution
  theta = parameters[sample,"GenomeCov_NB_theta"]
  mu = parameters[sample,"GenomeCov_NB_mu"]
  genome_expected_coverage = rnbinom (n = n_test, size = theta, mu = mu)
  plot(density(genome_simulated_coverage), col="red", 
       lwd=2, main = paste(sample," expected vs observed coverage distribution"))
  lines (density(genome_expected_coverage), col="blue", lwd = 2, add = TRUE)
  legend(1, 95, legend=c("Simulation", "Real data"),
         col=c("red", "blue"), cex=0.8)
  # remove tests with NA
  FisherTest = as.numeric(as.character(FisherTest[!is.na(FisherTest[,1]),1]))
  hist(FisherTest, breaks=100, freq = F)
  # Apply FDR
  FisherAdj = p.adjust(FisherTest, method = "fdr")
  
  # output metrics
  cutoff = 0.05
  percentage_FET_positives = 100 * length(FisherAdj[FisherAdj<cutoff]) / length(FisherAdj)
  count_FET_positives = length(FisherAdj[FisherAdj<cutoff])
  n_loci_tested = length(FisherAdj)
  
  # Results formatting
  res = c(round(parameters[sample,"GenomeCov_NB_mu"],0), 
          cutoff, 
          percentage_FET_positives,
          n_loci_tested, 
          count_FET_positives,oithona[sample],
          oithona[sample]-percentage_FET_positives,
          100*(oithona[sample]-percentage_FET_positives)/oithona[sample]
  )
  results = rbind(results, res)
}

colnames (results) = c( "Genome coverage", "q-value cutoff", 
                        "%FET positives","Number of tested loci",
                        "counts","%FET real data",
                        "% true psADE in real data", 
                        "% of detected psADE - noise")

rownames(results) = rownames(parameters)
write.table(results,"simulation_FET_results.txt", sep = "\t", quote=FALSE)





