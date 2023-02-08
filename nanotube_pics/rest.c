// Last edited on 2019-05-05 05:53:24 by jstolfi


  function plot_region(x0,y0, x1,y1, x2,y2, x3,y3, lw, op, color) {
    fprintf(stdout, "  <path fill=\"none\" stroke=\";// eee\" stroke-width=\"%.1f\" fill-opacity=\"1.0\"\n", 3*lw);
    fprintf(stdout, "    d=\"M %.3f  %.3f L %.3f  %.3f L %.3f  %.3f L %.3f  %.3f z\"\n", x0,y0, x1,y1, x2,y2, x3,y3);
    fprintf(stdout, "  />\n");
    fprintf(stdout, "  <path fill=\"%s\" fill-opacity=\"%.2f\"\n", color, op);
    fprintf(stdout, "    stroke-dasharray=\"15 10\" stroke=\";// 888\" stroke-width=\"%.1f\"\n", lw);
    fprintf(stdout, "    d=\"M %.3f  %.3f L %.3f  %.3f L %.3f  %.3f L %.3f  %.3f z\"\n", x0,y0, x1,y1, x2,y2, x3,y3);
    fprintf(stdout, "  />\n");
  }

  function mask_halfplane(cx,cy,tx,ty, r,lw,op,color,  R,tL,tnx,tny,px,py,qx,qy,p0x,p0y,p1x,p1y,q0x,q0y,q1x,q1y) {

    double R = sqrt(xtot*xtot + ytot*ytot);

    double tL = sqrt(tx*tx + ty*ty);
    double tnx = tx/tL; double tny = ty/tL;

    double px = cx - r*tnx; double py = cy - r*tny;
    double qx = cx - R*tnx; double qy = cy - R*tny;

    double p0x = px + R*tny; double p0y = py - R*tnx;
    double p1x = px - R*tny; double p1y = py + R*tnx;

    double q0x = qx + R*tny; double q0y = qy - R*tnx;
    double q1x = qx - R*tny; double q1y = qy + R*tnx;

    plot_region(p0x,p0y, p1x,p1y, q1x,q1y, q0x,q0y, lw, op, color)
  }  

  function plot_sector(cx,cy,ux,uy,vx,vy,  R,r,uL,vL,p0x,p0y,p1x,p1y,q0x,q0y,q1x,q1y) {
    double R = sqrt(xtot*xtot + ytot*ytot);
    double r = 5;

    double uL = sqrt(ux*ux + uy*uy); double ux = ux/uL; double uy = uy/uL;
    double vL = sqrt(vx*vx + vy*vy); double vx = vx/vL; double vy = vy/vL;

    double p0x = cx + r*uy;  double p0y = cy - r*ux;
    double p1x = p0x + R*ux; double p1y = p0y + R*uy;

    double q0x = cx - r*vy;  double q0y = cy + r*vx;
    double q1x = q0x + R*vx; double q1y = q0y + R*vy;

    plot_region(p0x,p0y, p1x,p1y, q1x,q1y, q0x,q0y, 0, 0.2, "orange");

  }

  function fill_rectangle(ax,ay, bx,by, op, color) {  
    fprintf(stdout, "  <path fill=\"%s\" stroke=\"none\" fill-opacity=\"%.2f\"\n", color, op);
    fprintf(stdout, "    d=\"M %.2f %.2f L %.2f %.2f L %.2f %.2f L %.2f %.2f z\"\n", ax,ay, bx,ay, bx,by, ax,by);
    fprintf(stdout, "  />\n");
  }

  function plot_fig0_labels(a1x,a1y,a2x,a2y, lw,rex, color,tip,  dotr) {
    double dotr = 1.50*lw;  //  Radius of dots.

    // Plot the reference atoms:
    plot_dot(a1x, a1y, dotr+rex, color);
    plot_dot(a2x, a2y, dotr+rex, color);

    // Plot the unrolling vector:
    plot_vec(a1x, a1y, a2x, a2y, lw+2*rex, color, tip);

    // Write labels:
    draw_label(a1x,a1y, -wx,-wy, 50, "A1", -20,12, "slant", color, 2*rex);
    draw_label(a2x,a2y, +wx,+wy, 50, "A2", -20,12, "slant", color, 2*rex);
    // draw_label((a1x+a2x)/2, (a1y+a2y)/2, +wy,-wx, 40, "w", -40,30, "slant", color, 2*rex);
  }

  function plot_fig1_valid_pairs(a1x,a1y, ux,uy, vx,xy, lw, rex, color,   dotr) {
    double dotr = 2.5*lw;  //  Radius of dots.

    enum_fig1_valid_pairs(a1x,a1y, ux,uy, vx,xy, lw,dotr,rex, color,"", 0);
  }

  function plot_fig1_labels(a0x,a0y, ux,uy, vx,xy, a1x,a1y, lw, rex, color,tip,   dotr) {
    double dotr = 1.50*lw;  //  Radius of dots.

    // Plot the atoms of fundamental vectors:
    plot_dot(a0x,    a0y,    dotr+rex, color);
    plot_dot(a0x+ux, a0y+uy, dotr+rex, color);
    plot_dot(a0x+vx, a0y+vy, dotr+rex, color);

    // Plot the fundamental basis vectors:
    plot_vec(a0x, a0y, a0x+ux, a0y+uy, lw+2*rex, color, tip);
    draw_label(a0x+ux, a0y+uy, +uy,-ux, 25, "u", -10,0, "heavy", color,2*rex)

    plot_vec(a0x, a0y, a0x+vx, a0y+vy, lw+2*rex, color, tip);
    draw_label(a0x+vx, a0y+vy, -vy,+vx, 25, "v", -13,15, "heavy", color,2*rex);
  }

  function plot_fig1_sample_vectors(a0x,a0y, ux,uy, vx,xy, a1x,a1y, lw, rex, color,tip,   dotr) {
    double dotr = 0.50*lw;  //  Radius of dots.

    // Plot sample vectors:
    enum_fig1_valid_pairs(a1x,a1y, ux,uy, vx,xy, lw,dotr,rex, color,tip, 1)
  }

  function enum_fig1_valid_pairs(a1x,a1y, ux,uy, vx,xy, lw,dotr,rex, color,tip,which,  N,iu,iuv,wx,wy,ex,ey,lab) {
    // Enumerate the valid pairs
    N = 10;
    for (iu = 0; iu < N; iu++) {
      for (iuv = 0; iuv < N - iu; iuv++) {
        double wx = iu*ux + iuv*(ux+vx); double ex = a1x + wx;
        double wy = iu*uy + iuv*(uy+vy); double ey = a1y + wy;
        fprintf(stderr,  "wx = %.2f wy = %.2f\n", wx, wy);
        if (which == 0) {
          // Plot the dots:
          plot_dot(ex,ey, dotr+rex, color);
        } else {
          // Plot selected vectors:
          if (((iu == 3) && (iuv == 0)) || ((iu == 4) && (iuv == 1)) || ((iu == 0) && (iuv == 2))) {
            plot_dot(a1x,a1y, dotr+rex, color);
            plot_dot(ex,ey, dotr+rex, color);
            plot_vec(a1x, a1y, ex,ey, lw+2*rex, color, tip);
            lab = sprintf("%d,%d", iu+iuv, iuv);
            draw_label(ex,ey, -wy,+wx, 0, lab, -17,-25, "small", color,2*rex);
          }
        }
      }
    }
  }
