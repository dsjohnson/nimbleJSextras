######## Figure 5 #######

tbeta = function(x, alpha1, beta1, beta2){
  
  cnst = beta(alpha1, beta1-beta2)/beta(alpha1, beta1)
  cnst*beta2*(1-x)^(beta2-1)*pbeta(x, alpha1, beta1-beta2)
  
}


rbeta_joint = function(u1, u2, a1, a2, b1, b2){
  
  cnst = (gamma(a1 + b1)*gamma(a2 + b2))/(gamma(a1)*gamma(b1)*gamma(a2)*gamma(b2))
  p1 = (u1^(a1-1))*((1-u1)^(b1-b2-a2))*((u2-u1)^(a2-1))*((1-u2)^(b2-1))
  cnst*p1*(u2>u1)*(u2<1)
  
}

u1grid = seq(0, 1, length.out = 250)
u2grid = seq(0, 1, length.out = 250)
pi_u1u2 = outer(u1grid, u2grid, rbeta_joint, a1 = 1.4, b1 = 5.6, a2 = 0.826, b2 = 2.478)
marg_u2 = colSums(pi_u1u2, na.rm = T)
marg_u1 = rowSums(pi_u1u2, na.rm = T)
marg_u2_norm = marg_u2/sum(marg_u2*(u2grid[2]-u2grid[1]))
marg_u1_norm = marg_u1/sum(marg_u1*(u1grid[2]-u1grid[1]))

######## FIGURE 5a #######

jpeg("Figure5a.jpg", width = 1200, height = 800, res = 150)

ggplot() +
  stat_function(fun = dbeta, args = list(shape1 = 0.5, shape2 = 1.5), color = I("firebrick"), linewidth = 1) +
  stat_function(fun = tbeta, args = list(alpha1 = 0.5, beta1 = 1.5, beta2 = .5), 
                color = I("darkorange"), linewidth = 1) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "u", y = "density") +
  theme_bw() +
  theme(text = element_text(size = 16))

dev.off()

######## FIGURE 5b #######

jpeg("Figure5b.jpg", width = 1200, height = 800, res = 150)

ggplot() +
  geom_line(aes(x = u1grid, y = marg_u1_norm), color = I("firebrick"), linewidth = 1) +
  geom_line(aes(x = u2grid, y = marg_u2_norm), color = I("darkorange"), linewidth = 1) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "u", y = "density") +
  theme_bw() +
  theme(text = element_text(size = 16))

dev.off()
