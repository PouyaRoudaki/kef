probs_test <- get_dens_or_prob(centered_kernel_mat_at_sampled, centered_kernel_mat_at_grid,
                          centered_kernel_self_grid,
                          sampled_x,
                          x_grid,
                          lambda_hat,
                          weights_hat,
                          type_of_p_is_prob = TRUE,
                          type_of_q_is_prob = TRUE,
                          method_of_p_calculation = "neighborhood_grid")

df1 <- data.frame(x1= sampled_x, y1 = colMeans(centered_kernel_mat_at_sampled))
df2 <- data.frame(x2= sampled_x, y2 = as.vector(probs_test$grid_x %*% t(centered_kernel_mat_at_grid)))
ggplot(df1,aes(x = x1,y=y1))+
  geom_point(colour = "red")+
  labs(x = "sampled_x",y = "Kernel Means")+
  theme_bw()

ggplot(df2,aes(x = x2,y=y2))+
  geom_point(colour = "blue")+
  labs(x = "sampled_x",y = "Kernel Mean Embedding")+
  theme_bw()

ggplot(df1,aes(x = x1,y=y1))+
  geom_point(colour = "red")+
  geom_point(data = df2,aes(x = x2,y=y2), colour = "blue")+
  labs(x = "sampled_x",y = "The Difference")+
  theme_bw()






summary(centered_kernel_mat_at_sampled)
