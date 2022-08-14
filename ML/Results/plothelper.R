library(ggplot2)
library(wesanderson)
library(plotly)
library(dplyr)


#figure 1A
plotly::ggplotly(
  ggplot(data = df, aes(x = accuracy, y=classifier,
                        fill=type))+
    geom_bar(stat="identity", color="black", position=position_dodge())+
    facet_grid(rows = vars(age), cols = vars(genes),
               scales="free_y", switch = 'y')+
    scale_fill_manual(values=wes_palette(n=2, name="Darjeeling2"),
                      name="Data Set", labels=c("Testing", "Training"))+
    theme_classic() +
    labs(title="Performance of Multiple Classifiers after RFE",
         x="",
         y="") +
    theme(plot.title=element_text(size=22, face="bold",hjust = 0.5, vjust=0.2,
                                  margin = margin(t = 0, r = 10, b = 10, l = 0)),
          axis.text.x=element_text(size=14, face="bold"),
          axis.text.y=element_text(size=14, face="bold"),
          axis.title.x=element_text(size=18, face="bold",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y=element_text(size=18, face="bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text=element_text(size=20, face="bold"),
          legend.title =element_text(size=18, face="bold"),
          strip.text.x = element_text(size=18, colour = "red"),
          strip.text.y = element_text(size=18, colour = "red"))+
    scale_y_discrete(limits=rev)+
    scale_x_continuous(labels = c(0,0.25,0.50,0.75,1))
) %>%
  layout(legend = list(orientation = "h", x =0.3, y =-0.12),
         title = list(y=0.98, x=0.15, xref="plot"),
         margin = list(l = 170, t = 100, r= 50),
         xaxis = list(title = list(text ='',
                                   standoff=20)),
         yaxis = list(title = list(text ='',
                                   standoff=30)),
         plot_bgcolor='white')


#figure 1B

plotly::ggplotly(
  ggplot(mRNA, aes(x=features)) +
    geom_line(aes(y = Train_accuracy, colour = "Training Accuracy"),
              size=1, linetype = "solid") +
    geom_point(aes(y = Train_accuracy, colour = "Training Accuracy"),
               size=2)+
    geom_line(aes(y = Test_accuracy, colour = "Testing Accuracy"),
              size=1, linetype = "solid") +
    geom_point(aes(y = Test_accuracy, colour = "Testing Accuracy"),
               size=2)+
    geom_line(aes(y = precision, colour = "Precision Score"),
              size=1, linetype = "solid") +
    geom_point(aes(y = precision, colour = "Precision Score"),
               size=2)+
    geom_line(aes(y = recall, colour = "Recall Score"),
              size=1, linetype = "solid") +
    geom_point(aes(y = recall, colour = "Recall Score"),
               size=2)+
    geom_line(aes(y = f1, colour = "F1 Score"),
              size=1, linetype = "solid") +
    geom_point(aes(y = f1, colour = "F1 Score"),
               size=2)+
    facet_grid(rows = vars(Age), cols = vars(classifier),
               scales="free", switch = 'both')+
    theme_bw() +
    labs(title="Model performance after hyperparameter tuning - mRNA dataset",
         x="Number of Features",
         y="Score",
         colour="Model Scoring Metrics") +
    theme(plot.title=element_text(size=20, face="bold",hjust = 0.6,
                                  margin = margin(t = 0, r = 10, b = 10, l = 0)),
          axis.text.x=element_text(size=16, face="bold"),
          axis.text.y=element_text(size=16, face="bold"),
          axis.title.x=element_text(size=20, face="bold",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y=element_text(size=20, face="bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text=element_text(size=16, face="bold"),
          legend.title =element_text(size=16, face="bold"),
          strip.text.x = element_text(size=16, colour = "red"),
          strip.text.y = element_text(size=16, colour = "red"))+
            scale_x_continuous(breaks= seq(0,90,by=15), 
                               limits = c(0,95), expand = c(0,0))
) %>% layout(legend = list(orientation = "h", x =0, y =-0.15),
         title = list(y=0.96, x=0.13, xref="plot"),
         margin = list(l = 100, t = 100, r= 50, b=100),
         plot_bgcolor='White')




plotly::ggplotly(
  ggplot(miRNA, aes(x=features)) +
    geom_line(aes(y = Train_accuracy, colour = "Training Accuracy"),
              size=1, linetype = "solid") +
    geom_point(aes(y = Train_accuracy, colour = "Training Accuracy"),
               size=2)+
    geom_line(aes(y = Test_accuracy, colour = "Testing Accuracy"),
              size=1, linetype = "solid") +
    geom_point(aes(y = Test_accuracy, colour = "Testing Accuracy"),
               size=2)+
    geom_line(aes(y = precision, colour = "Precision Score"),
              size=1, linetype = "solid") +
    geom_point(aes(y = precision, colour = "Precision Score"),
               size=2)+
    geom_line(aes(y = recall, colour = "Recall Score"),
              size=1, linetype = "solid") +
    geom_point(aes(y = recall, colour = "Recall Score"),
               size=2)+
    geom_line(aes(y = f1, colour = "F1 Score"),
              size=1, linetype = "solid") +
    geom_point(aes(y = f1, colour = "F1 Score"),
               size=2)+
    facet_grid(rows = vars(Age), cols = vars(classifier),
               scales="free", switch = 'both')+
    theme_bw() +
    labs(title="Model performance after hyperparameter tuning - miRNA dataset",
         x="Number of Features",
         y="Score",
         colour="Model Scoring Metrics") +
    theme(plot.title=element_text(size=20, face="bold",hjust = 0.6,
                                  margin = margin(t = 0, r = 10, b = 10, l = 0)),
          axis.text.x=element_text(size=16, face="bold"),
          axis.text.y=element_text(size=16, face="bold"),
          axis.title.x=element_text(size=20, face="bold",
                                    margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y=element_text(size=20, face="bold",
                                    margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text=element_text(size=16, face="bold"),
          legend.title =element_text(size=16, face="bold"),
          strip.text.x = element_text(size=16, colour = "red"),
          strip.text.y = element_text(size=16, colour = "red"))+
    scale_x_continuous(breaks= seq(0,100,by=20), 
                       limits = c(0,110), expand = c(0,0))
)%>%
  layout(legend = list(orientation = "h", x =0, y =-0.15),
         title = list(y=0.96, x=0.13, xref="plot"),
         margin = list(l = 100, t = 100, r= 50, b=100),
         plot_bgcolor='White')

