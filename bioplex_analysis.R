setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/BioRad BioPlex Assay/")
library(car)
library(ggpubr)
library(dunn.test)
library(ggpubr)
library(ggplot2)
rm(list=ls())

data1 = read.csv("RvE1 JAX Luminex Diabetes Panel.csv", header = TRUE, fill = TRUE, row.names = 1)
columns <- colnames(data1[,2:10]) #take all columns except for the diet group column

sink('Luminex_RvE1_JAX_statistics_2.19.20.txt') #open a text file to output all the stats in

for(col in columns){ #only loop through the column names specified above
  print(col) #printing column name
  formula = as.formula(paste(col, "~ Groups")) #creating the comparison formula with the column
  
  if(shapiro.test(data1[,col])$p.value > 0.05){ #if data is normal
    if(leveneTest(formula, data1)$Pr[[1]] > 0.05){ #if the data has equal variances
      my_anova <- aov(formula, data= data1)
      Anova(my_anova, type = "III") #must run a type III ANOVA because the sample size is not balanced
      print("Type 3 One-Way ANOVA")
      print(TukeyHSD(my_anova)) #Tukey adjusts the p-values (TukeyHSD can only be done after an ANOVA)
    }else { #if the data is normal but does NOT have equal variances
      #run the Welch ANOVA
      print("One-Way Welch ANOVA")
      print(oneway.test(formula, data = data1))
      #run a Pairwise t-test with no assumption of equal variances
      print("Pairwise T-test (no assumption of equal variance)")
      print(pairwise.t.test(data1[,col], data1$Groups, pool.sd = FALSE))
    }
  } else { #if the data is not normal
    if(leveneTest(formula, data1)$Pr[[1]] > 0.05){ #data is equal variance
      #dunns/kruskal wallis code
      print("Dunn's Test")
      print(dunn.test(data1[,col], data1$Groups))
    } else{ #data is not equal variance
      print("Wilcoxon Test (non-paired)")
      print(compare_means(formula,  data = data1, method = "wilcox.test", paired = FALSE))
    }
  }
}

#output results to a text file
sink()


#Make multiple plots (no stats annotation)

pdf("Luminex_RvE1_JAX_plots_test.pdf")

for(col in columns){
  print(ggbarplot(data1, x = "Groups", y = col, 
            add = c("mean_se", "point"), order = c("Con", "HF+Veh", "HF+RvE1"),
            color = "Groups", palette = "aaas", title = col,
            position = position_dodge(0.8)))
}

dev.off()

#############custom plots##########
ggbarplot(data1, x = "Groups", y = "Insulin", 
          add = c("mean_se", "point"), order = c("Con", "HF+Veh", "HF+RvE1"),
          color = "Groups", palette = "aaas", title = "Insulin",
          position = position_dodge(0.8))

ggbarplot(data1, x = "Groups", y = "Leptin", 
          add = c("mean_se", "point"), order = c("Con", "HF+Veh", "HF+RvE1"),
          color = "Groups", palette = "aaas", title = "Leptin",
          position = position_dodge(0.8))

ggbarplot(data1, x = "Groups", y = "Glucagon", 
          add = c("mean_se", "point"),order = c("Con", "HF+Veh", "HF+RvE1"),
          color = "Groups", palette = "aaas", title = "Glucagon",
          position = position_dodge(0.8))

ggbarplot(data1, x = "Groups", y = "Ghrelin", 
          add = c("mean_se", "point"),order = c("Con", "HF+Veh", "HF+RvE1"),
          color = "Groups", palette = "aaas", title = "Ghrelin",
          position = position_dodge(0.8))

ggbarplot(data1, x = "Groups", y = "Adiponectin", 
          add = c("mean_se", "point"),order = c("Con", "HF+Veh", "HF+RvE1"),
          color = "Groups", palette = "aaas", title = "Adiponectin",
          position = position_dodge(0.8))

ggbarplot(data1, x = "Groups", y = "GIP", 
          add = c("mean_se", "point"),order = c("Con", "HF+Veh", "HF+RvE1"),
          color = "Groups", palette = "aaas", title = "GIP",
          position = position_dodge(0.8))

ggbarplot(data1, x = "Groups", y = "GLP.1", 
          add = c("mean_se", "point"),order = c("Con", "HF+Veh", "HF+RvE1"),
          color = "Groups", palette = "aaas", title = "GLP.1",
          position = position_dodge(0.8))

ggbarplot(data1, x = "Groups", y = "PAI.1", 
          add = c("mean_se", "point"),order = c("Con", "HF+Veh", "HF+RvE1"),
          color = "Groups", palette = "aaas", title = "PAI.1",
          position = position_dodge(0.8))

ggbarplot(data1, x = "Groups", y = "Resistin", 
          add = c("mean_se", "point"),order = c("Con", "HF+Veh", "HF+RvE1"),
          color = "Groups", palette = "aaas", title = "Resistin",
          position = position_dodge(0.8))
#############
