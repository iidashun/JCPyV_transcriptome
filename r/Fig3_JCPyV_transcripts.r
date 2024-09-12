# 2023.1.20 originally created by Sotaro Mine
# 2023.12.13 modified by Shun Iida

# This program visualizes following transcripts.
# 1. Annotated transcripts
# 2. Non-annotated transcripts; >= 0.1% of early/late strands in both IMR-32 and 293
# 3. Non-annotated transcripts; wrapatound transcripts shown in Figs 4 and 5 (including transcripts < 0.1% of early/late strands)


##### settings #####

# load library
library(circlize)

# open tiff device
tiff("outputs/figs/fig3/JCPyV_transcripts.tiff", width = 5000, height = 5000, units = 'px', res = 300)

# change parameters as follows
circos.par(start.degree = 90, clock.wise = FALSE, xaxis.clock.wise = FALSE)

# define functions
# to indicate splicing acceptors/donors and polyA site
draw_radial_line_and_text <- function(df_splice) {
  for (i in 1:nrow(df_splice)) {
    bp <- df_splice$bp[i]
    type <- df_splice$type[i] # d: donor, a: acceptor, w: dual donor/acceptor, p: polyA
    bp.degree <- ((360 * bp) / 5130) + 90
    arrow_color <- switch(
      type,
      d = "#87CEEB",
      a = "#BEBEBE",
      w = "#B42DB4",
      p = "#000000"
    )
  
    if(bp %in% c(400, 492)) {
      arrows(0, 0, (0.507 * cos(bp.degree / 180 * pi)), (0.507 * sin(bp.degree / 180 * pi)), length = 0.001, lwd = 1.5, col = arrow_color)
      arrows((0.695 * cos(bp.degree / 180 * pi)), (0.695 * sin(bp.degree / 180 * pi)), cos(bp.degree / 180 * pi), sin(bp.degree / 180 * pi), length = 0.001, lwd = 1.5, col = arrow_color)
    } else {
      arrows(0, 0, cos(bp.degree / 180 * pi), sin(bp.degree / 180 * pi), length = 0.001, lwd = 1.5, col = arrow_color)
    }

    text(ifelse(type == "p", 1.025 * cos(((bp.degree) * pi)/ 180), ifelse(bp < 2561, 0.975 * cos(((bp.degree - 0.9) * pi)/ 180), ifelse(bp > 2561, 0.975 * cos(((bp.degree + 0.9) * pi)/ 180), 0))),
         ifelse(type == "p", 1.025 * sin(((bp.degree) * pi)/ 180), ifelse(bp < 2561, 0.975 * sin(((bp.degree - 0.9) * pi)/ 180), ifelse(bp > 2561, 0.975 * sin(((bp.degree + 0.9) * pi)/ 180), 0))),
         ifelse(type == "p", "Poly(A)", bp),
         cex = 1.4,
         col = "#000000",
         srt = ifelse(type == "p", 0, ifelse(bp < 2561, bp.degree + 180, ifelse(bp > 2561, bp.degree, 0)))
        )
  }
}

# to indicate JCPyV ORFs
draw_JCPyV_ORF <- function(start, end, label, color, textcolor, track) {
  circos.rect(start, ((track - 1) / 3), end, ((track) / 3), col = color)
  circos.text((start + end) / 2, (((track) /3) + ((track -1) /3)) /2, label, facing = "bending.inside", cex = 1.0, col = textcolor)
}

# initializing; JCPyV genome = 5,130 bp
circos.initialize("JCPyV", xlim = c(0, 5130), ring = TRUE)

# plotting from outer tracks to innere tracks

df = data.frame(
  name = c(" "), start = c(0), end = c(5130))


##### plotting transcripts #####

# non-annotated transcripts
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.5)

#################### 1st ####################

# IMR-L2, IMR-L10, IMR-L92; wraparoud
# If 298-400 is spliced out; IMR-L32, IMR-L67
circos.lines(c(267, 298), c(0.00, 0.00), lwd = 1.5, col = "#FF0000")
circos.lines(c(298, 400), c(0.00, 0.00), lwd = 1.5, col = "#FF0000")
circos.lines(c(400, 492), c(0.01, 0.01), lwd = 1.5, col = "#FF0000")
circos.lines(c(400, 492), c(-0.01, -0.01), lwd = 1.5, col = "#FF0000")
circos.lines(c(492, 1426), c(0.00, 0.00), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(1426, 2556), c(0.00, 0.00), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.00, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.0, "VP1-401", facing = "bending.inside", cex = 1.3, col = "#000000")
circos.text((400 + 492) / 2, 0.05, "(x2-x4)", facing = "bending.inside", cex = 1.3, col = "#000000")

# IMR-E6
circos.lines(c(4273, 5023), c(0.00, 0.00), lwd = 1.5, col = "#FF0000")
circos.lines(c(2918, 4273), c(0.00, 0.00), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(2566, 2918), c(0.00, 0.00), lwd = 1.5, col = "#FF0000")
circos.arrow(2566, 2566, y = 0.00, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), arrow.position = "start", border = NA, col = "#FF0000")
circos.text(5115, 0.0, "E6", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 2nd ####################

# IMR-L3, IMR-L14; wraparoud
# If 298-376 is spliced out; 293-L60 
circos.lines(c(267, 376), c(0.10, 0.10), lwd = 1.5, col = "#FF0000")
circos.lines(c(376, 492), c(0.11, 0.11), lwd = 1.5, col = "#FF0000")
circos.lines(c(376, 492), c(0.09, 0.09), lwd = 1.5, col = "#FF0000")
circos.lines(c(492, 1426), c(0.10, 0.10), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(1426, 2556), c(0.10, 0.10), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.10, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.1, "VP1-377", facing = "bending.inside", cex = 1.3, col = "#000000")
circos.text((400 + 492) / 2, 0.15, "(x2-x3)", facing = "bending.inside", cex = 1.3, col = "#000000")

# IMR-E9
circos.lines(c(4493, 5023), c(0.10, 0.10), lwd = 1.5, col = "#FF0000")
circos.lines(c(4426, 4493), c(0.10, 0.10), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(4273, 4426), c(0.10, 0.10), lwd = 1.5, col = "#FF0000")
circos.lines(c(2777, 4273), c(0.10, 0.10), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(2566, 2777), c(0.10, 0.10), lwd = 1.5, col = "#FF0000")
circos.arrow(2566, 2566, y = 0.10, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), arrow.position = "start", border = NA, col = "#FF0000")
circos.text(5115, 0.1, "E9", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 3rd ####################

# IMR-L5, IMR-L22, IMR-L46; wraparoud
# If 298-400 is spliced out; IMR-L48
circos.lines(c(267, 400), c(0.20, 0.20), lwd = 1.5, col = "#FF0000")
circos.lines(c(400, 492), c(0.21, 0.21), lwd = 1.5, col = "#FF0000")
circos.lines(c(400, 492), c(0.19, 0.19), lwd = 1.5, col = "#FF0000")
circos.lines(c(492, 2556), c(0.20, 0.20), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.20, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.2, "VP2/3-401", facing = "bending.inside", cex = 1.3, col = "#000000")
circos.text((400 + 492) / 2, 0.25, "(x2-x4)", facing = "bending.inside", cex = 1.3, col = "#000000")

# IMR-E11
circos.lines(c(4493, 5023), c(0.20, 0.20), lwd = 1.5, col = "#FF0000")
circos.lines(c(4426, 4493), c(0.20, 0.20), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(4273, 4426), c(0.20, 0.20), lwd = 1.5, col = "#FF0000")
circos.lines(c(2704, 4273), c(0.20, 0.20), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(2566, 2704), c(0.20, 0.20), lwd = 1.5, col = "#FF0000")
circos.arrow(2566, 2566, y = 0.20, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), arrow.position = "start", border = NA, col = "#FF0000")
circos.text(5115, 0.2, "E11", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 4th ####################

# IMR-L7, IMR-L37; wraparoud
# If 298-376 is spliced out; IMR-L102
circos.lines(c(267, 376), c(0.30, 0.30), lwd = 1.5, col = "#FF0000")
circos.lines(c(376, 492), c(0.31, 0.31), lwd = 1.5, col = "#FF0000")
circos.lines(c(376, 492), c(0.29, 0.29), lwd = 1.5, col = "#FF0000")
circos.lines(c(492, 2556), c(0.30, 0.30), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.30, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.3, "VP2/3-377", facing = "bending.inside", cex = 1.3, col = "#000000")
circos.text((400 + 492) / 2, 0.35, "(x2-x3)", facing = "bending.inside", cex = 1.3, col = "#000000")

# IMR-E12
circos.lines(c(4273, 5023), c(0.30, 0.30), lwd = 1.5, col = "#FF0000")
circos.lines(c(2777, 4273), c(0.30, 0.30), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(2566, 2777), c(0.30, 0.30), lwd = 1.5, col = "#FF0000")
circos.arrow(2566, 2566, y = 0.30, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), arrow.position = "start", border = NA, col = "#FF0000")
circos.text(5115, 0.3, "E12", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 5th ####################

# IMR-L9
circos.lines(c(267, 298), c(0.40, 0.40), lwd = 1.5, col = "#FF0000")
circos.lines(c(298, 400), c(0.40, 0.40), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(400, 492), c(0.40, 0.40), lwd = 1.5, col = "#FF0000")
circos.lines(c(492, 1426), c(0.40, 0.40), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(1426, 2556), c(0.40, 0.40), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.40, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.4, "L9", facing = "bending.inside", cex = 1.3, col = "#000000")

# IMR-E13
circos.lines(c(4493, 5023), c(0.40, 0.40), lwd = 1.5, col = "#FF0000")
circos.lines(c(4426, 4493), c(0.40, 0.40), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(4273, 4426), c(0.40, 0.40), lwd = 1.5, col = "#FF0000")
circos.lines(c(2918, 4273), c(0.40, 0.40), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(2566, 2918), c(0.40, 0.40), lwd = 1.5, col = "#FF0000")
circos.arrow(2566, 2566, y = 0.40, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), arrow.position = "start", border = NA, col = "#FF0000")
circos.text(5115, 0.4, "E13", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 6th ####################

# IMR-L13
circos.lines(c(267, 298), c(0.50, 0.50), lwd = 1.5, col = "#FF0000")
circos.lines(c(298, 400), c(0.50, 0.50), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(400, 2556), c(0.50, 0.50), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.50, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.5, "L13", facing = "bending.inside", cex = 1.3, col = "#000000")

# 293-E16, IMR-E15, IMR-E14, IMR-E16; wraparoud
circos.lines(c(4426, 5023), c(0.50, 0.50), lwd = 1.5, col = "#FF0000")
circos.lines(c(4273, 4426), c(0.51, 0.51), lwd = 1.5, col = "#FF0000")
circos.lines(c(4273, 4426), c(0.49, 0.49), lwd = 1.5, col = "#FF0000")
circos.lines(c(2566, 4273), c(0.50, 0.50), lwd = 1.5, col = "#FF0000")
circos.arrow(2566, 2566, y = 0.5, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), arrow.position = "start", border = NA, col = "#FF0000")
circos.text(5115, 0.5, "SuperT", facing = "bending.inside", cex = 1.3, col = "#000000")
circos.text((4273 + 4426) / 2, 0.55, "(x3-x9)", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 7th ####################

# IMR-L17
circos.lines(c(267, 298), c(0.60, 0.60), lwd = 1.5, col = "#FF0000")
circos.lines(c(298, 521), c(0.60, 0.60), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(521, 2556), c(0.60, 0.60), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.60, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.6, "L17", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### Novel ####################

# Novel
circos.text(3300, 0.55, "Novel", facing = "bending.outside", cex = 1.3, col = "#000000")

#################### 8th ####################

# IMR-L18
circos.lines(c(267, 315), c(0.70, 0.70), lwd = 1.5, col = "#FF0000")
circos.lines(c(315, 400), c(0.70, 0.70), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(400, 492), c(0.70, 0.70), lwd = 1.5, col = "#FF0000")
circos.lines(c(492, 1426), c(0.70, 0.70), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(1426, 2556), c(0.70, 0.70), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.70, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.7, "L18", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 9th ####################

# IMR-L21
circos.lines(c(267, 400), c(0.80, 0.80), lwd = 1.5, col = "#FF0000")
circos.lines(c(400, 521), c(0.80, 0.80), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(521, 2556), c(0.80, 0.80), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.80, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.8, "L21", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 10th ####################

# IMR-L27
circos.lines(c(267, 1585), c(0.90, 0.90), lwd = 1.5, col = "#FF0000")
circos.lines(c(1585, 1953), c(0.90, 0.90), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(1953, 2556), c(0.90, 0.90), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.90, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.9, "L27", facing = "bending.inside", cex = 1.3, col = "#000000")

set_track_gap(mm_h(0))


# annotated transcripts
circos.track(ylim = c(0, 0.6), bg.border = NA, bg.col = "#FFFFE5", track.height = 0.25)

#################### 1st ####################

# M2
circos.lines(c(267, 492), c(0.0, 0.0), lwd = 1.5, col = "#FF0000")
circos.lines(c(492, 1426), c(0.0, 0.0), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(1426, 2556), c(0.0, 0.0), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.0, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.0, "M2", facing = "bending.inside", cex = 1.3, col = "#000000")

# LT
circos.lines(c(4770, 5023), c(0.0, 0.0), lwd = 1.5, col = "#FF0000")
circos.lines(c(4426, 4770), c(0.0, 0.0), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(2566, 4426), c(0.0, 0.0), lwd = 1.5, col = "#FF0000")
circos.arrow(2566, 2566, y = 0.0, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), arrow.position = "start", border = NA, col = "#FF0000")
circos.text(5115, 0.0, "LT", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 2nd ####################

# M3
circos.lines(c(267, 492), c(0.1, 0.1), lwd = 1.5, col = "#FF0000")
circos.lines(c(492, 1426), c(0.1, 0.1), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(1426, 1585), c(0.1, 0.1), lwd = 1.5, col = "#FF0000")
circos.lines(c(1585, 1953), c(0.1, 0.1), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(1953, 2556), c(0.1, 0.1), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.1, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.1, "M3", facing = "bending.inside", cex = 1.3, col = "#000000")

# ST
circos.lines(c(4493, 5023), c(0.1, 0.1), lwd = 1.5, col = "#FF0000")
circos.lines(c(4426, 4493), c(0.1, 0.1), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(2566, 4426), c(0.1, 0.1), lwd = 1.5, col = "#FF0000")
circos.arrow(2566, 2566, y = 0.1, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), arrow.position = "start", border = NA, col = "#FF0000")
circos.text(5115, 0.1, "ST", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 3rd ####################

# M4
circos.lines(c(267, 492), c(0.2, 0.2), lwd = 1.5, col = "#FF0000")
circos.lines(c(492, 1953), c(0.2, 0.2), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(1953, 2556), c(0.2, 0.2), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.2, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.2, "M4", facing = "bending.inside", cex = 1.3, col = "#000000")

# T'135
circos.lines(c(4770, 5023), c(0.2, 0.2), lwd = 1.5, col = "#FF0000")
circos.lines(c(4426, 4770), c(0.2, 0.2), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(4273, 4426), c(0.2, 0.2), lwd = 1.5, col = "#FF0000")
circos.lines(c(2918, 4273), c(0.2, 0.2), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(2566, 2918), c(0.2, 0.2), lwd = 1.5, col = "#FF0000")
circos.arrow(2566, 2566, y = 0.2, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), arrow.position = "start", border = NA, col = "#FF0000")
circos.text(5115, 0.2, "T'135", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 4th ####################

# M5
circos.lines(c(267, 732), c(0.3, 0.3), lwd = 1.5, col = "#FF0000")
circos.lines(c(732, 1426), c(0.3, 0.3), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(1426, 2556), c(0.3, 0.3), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.3, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.3, "M5", facing = "bending.inside", cex = 1.3, col = "#000000")

# T'136
circos.lines(c(4770, 5023), c(0.3, 0.3), lwd = 1.5, col = "#FF0000")
circos.lines(c(4426, 4770), c(0.3, 0.3), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(4273, 4426), c(0.3, 0.3), lwd = 1.5, col = "#FF0000")
circos.lines(c(2777, 4273), c(0.3, 0.3), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(2566, 2777), c(0.3, 0.3), lwd = 1.5, col = "#FF0000")
circos.arrow(2566, 2566, y = 0.3, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), arrow.position = "start", border = NA, col = "#FF0000")
circos.text(5115, 0.3, "T'136", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 5th ####################

# M6
circos.lines(c(267, 732), c(0.4, 0.4), lwd = 1.5, col = "#FF0000")
circos.lines(c(732, 1426), c(0.4, 0.4), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(1426, 1585), c(0.4, 0.4), lwd = 1.5, col = "#FF0000")
circos.lines(c(1585, 1953), c(0.4, 0.4), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(1953, 2556), c(0.4, 0.4), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.4, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.4, "M6", facing = "bending.inside", cex = 1.3, col = "#000000")

# T'165
circos.lines(c(4770, 5023), c(0.4, 0.4), lwd = 1.5, col = "#FF0000")
circos.lines(c(4426, 4770), c(0.4, 0.4), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(4273, 4426), c(0.4, 0.4), lwd = 1.5, col = "#FF0000")
circos.lines(c(2704, 4273), c(0.4, 0.4), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(2566, 2704), c(0.4, 0.4), lwd = 1.5, col = "#FF0000")
circos.arrow(2566, 2566, y = 0.4, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), arrow.position = "start", border = NA, col = "#FF0000")
circos.text(5115, 0.4, "T'165", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 6th ####################

# ORF1
circos.lines(c(267, 492), c(0.5, 0.5), lwd = 1.5, col = "#FF0000")
circos.lines(c(492, 521), c(0.5, 0.5), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(521, 1426), c(0.5, 0.5), lwd = 1.5, col = "#FF0000")
circos.lines(c(1426, 1585), c(0.51, 0.51), lwd = 1.5, col = "#FF0000")
circos.lines(c(1426, 1585), c(0.49, 0.49), lwd = 1.5, col = "#FF0000")
circos.lines(c(1585, 2556), c(0.5, 0.5), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.5, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.5, "ORF1", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### 7th ####################

# ORF2
circos.lines(c(267, 492), c(0.6, 0.6), lwd = 1.5, col = "#FF0000")
circos.lines(c(492, 521), c(0.6, 0.6), lwd = 1.5, col = "#0000FF", lty = "dashed")
circos.lines(c(521, 1426), c(0.6, 0.6), lwd = 1.5, col = "#FF0000")
circos.lines(c(1426, 1643), c(0.61, 0.61), lwd = 1.5, col = "#FF0000")
circos.lines(c(1426, 1643), c(0.59, 0.59), lwd = 1.5, col = "#FF0000")
circos.lines(c(1643, 2556), c(0.6, 0.6), lwd = 1.5, col = "#FF0000")
circos.arrow(2556, 2556, y = 0.6, width = 0.01, arrow.head.width = 0.05, arrow.head.length = cm_x(0.25), border = NA, col = "#FF0000")
circos.text(175, 0.6, "ORF2", facing = "bending.inside", cex = 1.3, col = "#000000")

#################### Annotated ####################

# Annotated
circos.text(3300, 0.55, "Annotated", facing = "bending.outside", cex = 1.3, col = "#000000")


# radial lines
splicing_site <- read.csv("inputs/JCPyV/csv/JCPyV_splicing_site.csv", header = TRUE, sep = ",")
draw_radial_line_and_text(splicing_site)

set_track_gap(mm_h(0))


##### plotting JCPyV ORFs #####

circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.165,
  panel.fun = function(x, y) {
    draw_JCPyV_ORF(883, 1560, "VP3", "#A5D5D8", "#000000", 1)
    draw_JCPyV_ORF(4495, 5013, "ST", "#F4777F", "#000000", 1)
    draw_JCPyV_ORF(526, 1560, "VP2", "#73A2C6", "#000000", 2)
    draw_JCPyV_ORF(1469, 2533, "VP1", "#00429D", "#FFFFFF", 3)
    draw_JCPyV_ORF(277, 492, "Agno", "#3CB371", "#000000", 3)
    draw_JCPyV_ORF(2603, 4426, "LT", "#93003A", "#FFFFFF", 3)
    draw_JCPyV_ORF(4771, 5013, "LT", "#93003A", "#FFFFFF", 3)
    circos.lines(c(0, 5130), c(1, 1), col = "#000000")
    circos.lines(c(4426, 4771), c(2/3, 2/3), col = "#000000", lty = "dashed")
  }
)

set_track_gap(mm_h(0))


##### virus name; JCPyV #####

draw.sector(center = c(0, 0), start.degree = 0, end.degree = 360, rou1 = 0.08, col = "#FFFFFF", border = "#FFFFFF")
text(0, 0, "JCPyV", cex = 1.3, font = 2)


##### device off #####

circos.clear() #これで設定をクリアする

dev.off()
