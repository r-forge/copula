## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


##------- FIXME{MM}: Shouldn't all these inherit from "htest" (and use its print() method)?
##        and ditto for all other ./*Tests.R  !!!

### EV test based on An CFG - An Pickands ######################################

evTestAA <- function(x, N = 1000,  derivatives = "Cn", m = 100)
{
  ## make pseudo-observations
  n <- nrow(x)
  u <- apply(x,2,rank)/(n+1)

  ## make grid
  g <- seq(1/m, 1 - 1/m, len = m)

  ## compute the test statistic
  s <- .C(evTestAA_stat,
          as.double(-log(u[,1])),
          as.double(-log(u[,2])),
          as.integer(n),
          as.double(g),
          as.integer(m),
          stat = double(1))$stat

  if (derivatives == "Cn")
    s0 <- .C(evTestAA,
             as.double(u[,1]),
             as.double(u[,2]),
             as.integer(n),
             as.double(g),
             as.integer(m),
             as.integer(N),
             s0 = double(N))$s0
  else
    s0 <- .C(evTestAA_derA,
             as.double(u[,1]),
             as.double(u[,2]),
             as.integer(n),
             as.double(g),
             as.integer(m),
             as.integer(N),
             s0 = double(N))$s0

  structure(class = "evTest",
            list(statistic=s, pvalue=(sum(s0 >= s)+0.5)/(N+1),s0=s0))
}

