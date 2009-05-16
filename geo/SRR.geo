// split-ring resonator structure
lc = 0.01;
w = 0.04;
a1 = 0.07;
a2 = 0.11;
b1 = 0.13;
b2 = 0.15;

Dx = 2*b2 + 0.03;
Nx = 6;
Dy = 0.15;
Ny = 4;
Dz = 2*b2 + 0.03;
Nz = 3;

Function Ring
  OrigNumber = newp;
  Point(OrigNumber) = {Origin[0], Origin[1], Origin[2], lc};
  Theta1 = Asin(w/2.0 / r1);
  pNumber[0] = newp;
  Point(pNumber[0]) = {Origin[0] + w/2.0, Origin[1], Origin[2] + sign * r1 * Cos(Theta1), lc};
  Theta2 = Asin(w/2.0 / r2);
  pNumber[1] = newp;
  Point(pNumber[1]) = {Origin[0] + w/2.0, Origin[1], Origin[2] + sign * r2 * Cos(Theta2), lc};
  pNumber[2] = newp;
  Point(pNumber[2]) = {Origin[0] + r1, Origin[1], Origin[2], lc};
  pNumber[3] = newp;
  Point(pNumber[3]) = {Origin[0] + r2, Origin[1], Origin[2], lc};
  pNumber[4] = newp;
  Point(pNumber[4]) = {Origin[0], Origin[1], Origin[2] - sign * r1, lc};
  pNumber[5] = newp;
  Point(pNumber[5]) = {Origin[0], Origin[1], Origin[2] - sign * r2, lc};
  // the lines
  elemNumber[0] = newreg;
  Line(elemNumber[0]) = {pNumber[0], pNumber[1]};
  elemNumber[1] = newreg;
  Line(elemNumber[1]) = {pNumber[2], pNumber[3]};
  elemNumber[2] = newreg;
  Line(elemNumber[2]) = {pNumber[4], pNumber[5]};
  // now the circles
  elemNumber[3] = newreg;
  Circle(elemNumber[3]) = {pNumber[0], OrigNumber, pNumber[2]};
  elemNumber[4] = newreg;
  Circle(elemNumber[4]) = {pNumber[1], OrigNumber, pNumber[3]};
  elemNumber[5] = newreg;
  Circle(elemNumber[5]) = {pNumber[2], OrigNumber, pNumber[4]};
  elemNumber[6] = newreg;
  Circle(elemNumber[6]) = {pNumber[3], OrigNumber, pNumber[5]};
  // now the surfaces
  loopNumber = newreg;
  Line Loop(loopNumber) = {elemNumber[3],elemNumber[1],-elemNumber[4],-elemNumber[0]};
  PlaneSurfNum[0] = newreg;
  Plane Surface(PlaneSurfNum[0]) = {loopNumber};
  loopNumber = newreg;
  Line Loop(loopNumber) = {elemNumber[1],elemNumber[6],-elemNumber[2],-elemNumber[5]};
  PlaneSurfNum[1] = newreg;
  Plane Surface(PlaneSurfNum[1]) = {loopNumber};
  // now the symmetries
  Symmetry {1,0,0,-Origin[0]} {
    Duplicata{ Surface{PlaneSurfNum[0], PlaneSurfNum[1]}; }
  }
Return

Function SRR
  r1 = a1; r2 = a2; sign = 1.0;
  Call Ring;
  r1 = b1; r2 = b2; sign = -1.0;
  Call Ring;
Return


For i In {0:Nx-1}
  For j In {0:Ny-1}
    For k In {0:Nz-1}
      Origin[] = {-(Nx-1) * Dx/2.0 + i*Dx, -(Ny-1) * Dy/2.0 + j*Dy, -(Nz-1) * Dz/2.0 + k*Dz};
      Call SRR;
    EndFor
  EndFor
EndFor
