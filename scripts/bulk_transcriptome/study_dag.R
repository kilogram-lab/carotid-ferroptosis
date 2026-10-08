# httpLOCAL_PROJECT_PATH


library("dagitty")
library(dagitty)
library(ggdag)
dag <- dagitty:: dagitty("dag{
                         X -> Z -> Y; X <- U -> Y;
                         {A R D J S} -> Y;
                         {B I} -> Z; 
                         I -> X; I -> Y;
                         D -> S; D -> J;
                         R -> Z;
                  }") 
