function visit_proc(e)

global saida;


o = Org(e);
d = Dest(e);
x = [o(1),d(1)];
y = [o(2),d(2)];

aresta.x = x;
aresta.y = y;

saida(length(saida)+1) = aresta

plot(x,y,'b-');