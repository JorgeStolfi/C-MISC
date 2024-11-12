/* 
  Last edited on 2009-02-10 08:30:43 by stolfi
  
  Programa para pre-processamento do dicionário de flexões de S. Carlos.
  
  UNDER DEVELOPMENT
  
  Criado por stolfi em 96-01-10, baseado no
  "trans.l" do tomasz de 95-06-08. 
  
  USO
  
    dfltovfl <entrada> <saida>
    
  onde <entrada> é o nome completo do dicionário de flexões de
  S. Carlos, no formato de quádruplas, e <saida> é o nome *sem extensão*
  dos arquivos de saída.  O programa grava dois arquivos:
    
    <saida>.vfa com o dicionário de flexões analítico, 
    <saida>.vfs com o dicionário de flexões sintético,
    <saida>.vms com as mensagens de erro.

  FUNÇÃO
 
  O dicionário de S. Carlos, depois de suficientemente 
  massageado, produz um conjunto de quuíntuplas ("entradas") 
  da forma
  
    (<palavra>,<ordem>,<classe>,<flexao>,<canonica>)
          
  onde  
    
    <canonica> é um string de até MAXLPAL caracteres do alfabeto,
        que é uma palavra da língua em sua flexão canônica.
  
    <classe> é uma letra que indica uma função gramatical
        geral que a <canonica> pode ter, como verbo, 
        substantivo, etc. (Vide tabela adiante)
      
    <flexao> é um string de até MAXNCAT-1 letras e dígitos que
        indicam uma flexão aplicável à <canonica> quando ela
        occorre como <classe>: gênero e número para substantivos,
        tempo e pessoa para verbos, etc.
       
    <palavra> é um string de até MAXLPAL caracteres do alfabeto,
        que representa a forma da <canonica> flexionada conforme
        indicado por <classe> e <flexao>. 
        
    <ordem> é um dígito de ordem que indica a importância
        relativa das várias entradas de uma mesma <palavra>.
 
  A concatenação de <classe> e <flexao> é chamada por nós de
  <categoria>.
       
  A <palavra> pode ser vazia, significando que essa particular
  <flexao> da forma <canonica> definidamente não existe. 
  
  O "vocabulário" é o conjunto de todas as <palavras> que aparecem
  neste arquivo.
  
 O dicionário de flexões será usado basicamente de três 
  maneiras: 
  
    busca "existencial": dada uma máscara,
    enumerar as <palavras> do vocabulário que satisfazem
    essa máscara. (Em particular, dada uma única <palavra>,
    verificar se ela está no vocabulário,)
  
    busca "sintética": dada uma tripla <canonica>, <classe>,
    <flexao>, determinar a <palavra> correspondente.
    
    busca "analítica": dada uma <palavra> do vocabulário, 
    determinar todas as quádruplas <ordem>,<canonica>,<classe>,<flexao>
    associadas a ela.  
    
  Para implementar essas buscas de maneira eficiente
  usando automatos, codificaremos o dicionário de flexões 
  por dois conjuntos de cadeias, representando aproximadamente a
  função flexionadora e sua inversa. Especificamente, para cada
  entrada que tem <palavra> não-vazia,  criamos
  duas cadeias ("sintética" e "analítica"):
  
    <palavra><a1><ordem><classe><flexao><a2><nsufpal><sufcanonica>#

    <canonica><s1><classe><flexao><s2><nsufcan><sufpalavra># 
          
  onde 

    <a1>, <a2>, <s1>, <s2> são delimitadores, 
        sendo que <s1> é distinto de <a1>;
     
    <sufcanonica> e <sufpalavra> são as terminações distintivas
        de <canonica> e <palavra>, respectivamente;
        
    <nsufpal> e <nsufcan> são os respectivos comprimentos.
   
  Mais precisamente, se <maxpref> é o maior prefixo comum entre 
  <palavra> e <canonica>, então <sufcanonica> e <sufpalavra> 
  são os respectivos sufixos, isto é,

    <canonica> = <maxpref><sufcanonica>
    <palavra>  = <maxpref><sufpalavra>
      
  Note que <sufcanonica> e/ou <sufpalavra> podem ser vazios.
  Os campos <nsufpal> <nsufcan> são dígitos hexatrisdecimais 
  (0-9a-z) que indicam o comprimento do outro sufixo:
   
    <nsufpal>  = |<sufpalavra>|
    <nsufcan>  = |<sufcanonica>|
      
  A versão corrente do programa supõe que estes comprimentos 
  são sempre menores ou iguais a 35.  
   
*/

#define _GNU_SOURCE
#include <stdio.h> 
#include <string.h>
#include <stdlib.h>

#define false 0
#define true 1 

typedef int bool_t;
typedef unsigned char byte;

/* Comprimentos maximos dos campos: */
#define MAXLPAL 34
#define MAXLORD  1
#define MAXLCLS 34
#define MAXLFLX 34
#define MAXLLIN 255
  
/* Max numero de erros: */
#define MAXERR 10000 
  
/* Max numero de registros a processar: */
#define MAXREG 1000000
        
/* Separadores de campos -- dicionário sintético de flexoes: */
#define SEPS1 ":"
#define SEPS2 "=" 

/* Separadores de campos -- dicionario analitico de flexoes: */
#define SEPA1 "="
#define SEPA2 ":"
  
long num_linhas = 0;        /* Conta linhas, incluindo vazias */
long num_vazias = 0;        /* Conta linhas vazias */
long num_erros = 0;         /* Conta linhas com formato invalido */
long num_inexistentes = 0;  /* Conta formas flexionadas inexistentes */
long num_palavras = 0;      /* Conta palavras do vocabulário */
long num_analiticas = 0;    /* Conta entradas do dic analítico de flexoes */
long num_sinteticas = 0;    /* Conta entradas do dic sintético de flexoes */

/* Arquivos */
static FILE 
  *dfl,   /* Dicionário de S. Carlos, formato quintuplas */
  *vfa,   /* Dicionário analítico de flexões */
  *vfs,   /* Dicionário sintético de flexões */
  *vms;   /* Mensagens */

#ifdef MSDOS
/* #include <process.h> */
#define RMODE "rt"
#define WMODE "wt"
#else
#define RMODE "r"
#define WMODE "w"
#endif

/* Internal prototypes: */

int main(int argc, char *argv[]);

FILE *openfile(char *name, char *mode);

bool_t le_linha(byte *lin, int maxlen, long *nlin); 
  /* 
    Lê proxima linha de "dfl", terminada por '\n'
    Devolve em "*nlin" o comprimento real da linha 
    (sem contar o '\n' final).  Coloca os bytes lidos
    em *lin (até o máximo de "maxlen" bytes, ignorando o excesso).
    Dá erro se a última linha não tiver '\n' no fim.
    Resultado é "false" sse "stdin" está esgotado.
  */
  
void processa_linha(byte *lin, long nlin);
  /* 
    Verifica, analisa, e processa uma linha,
    dados o comprimento total "nlin" e os primeiros MAXLLIN
    caracteres "lin[0..min(nlin,MAXLLIN)-1]".  Ambos excluem
    o NEWLINE final.
  */
  
bool_t testa_delim(byte **r, long *n, char del, byte *lin, long nlin);
  /*
    Verifica se o caracter "**r" é igual a "del"; em caso afirmativo
    incrementa "*r", decrementa "*n", e devolve "true";
    caso contrário imprime erro e devolve "false". 
  */
  
void mens_linha(char *msg, long num, byte *lin, long nlin); 
  /*
    Imprime uma mensagem "msg" referente à linha numero "num",
    cujo conteudo (truncado em MAXLLIN bytes) é "*lin" e 
    cujo comprimento é "*nlin".
  */
 
void mens_sintaxe(char *msg, long col, byte *lin, long nlin); 
  /*
    Imprime uma mensagem de erro de sintaxe "msg" referente 
    ao caracter indice "col" da linha numero "num_linhas",
    cujo conteudo (truncado em MAXLLIN bytes) é "*lin" e 
    cujo comprimento é "*nlin".
  */ 
  
void imprime_linha(byte *lin, long nlin);
  /* 
    Imprime a linha "lin" no arquivo "vms", com escapes e elipses,
    para fins de debugagem.
  */
 
void define_tabela(bool_t tb[], char *lets);
  /* 
    Define "tb[c]=true" para todo "c" em "*lets"; demais "false".
  */
  
bool_t cata(
    byte **buf, long *nbuf, 
    byte **res, long *nres,
    bool_t tb[], long maxres,
    byte *lin, long nlin
  );
  /* 
    Retorna em (*res) um apontador para o maior prefixo
    da cadeia (*rbuf)[0..(*nbuf)-1] que é composto de caracteres 
    da tabela "tb"; e devolve em *nres o comprimento desse prefixo.
    Também atualiza (*rbuf) e (*nbuf) para o restante da cadeia. 
    Se o prefixo tem comprimento maior  que "maxres", imprime
    uma mensagem de erro (com a linha "lin"), incrementa num_erros,
    e devolve "false"; caso contrário devolve "true".
  */
  
char *txtcat(char *a, char *b);
  /*
    Devolve uma nova cadeia que contém a concatenação de "a" 
    e "b".
  */
 
#define check(cond,msg) \
  { if(!(cond)) \
      { if (vms) fprintf(vms, "** %s\n", msg); \
        if (vms != stderr) fprintf(stderr, "** %s\n", msg); \
        exit(1); \
      }  \
  }

static bool_t tb_pal[256]; /* Bytes válidos em palavras */
static bool_t tb_cat[256]; /* Bytes válidos em dados gramaticais */ 
static bool_t tb_ord[256]; /* Bytes válidos em dígitos de ordem */ 

static char digit[36] = "0123456789abcdefghijklmnopqrstuvwxyz";

int main(int argc, char *argv[])
  { byte *lin = (byte*) malloc(MAXLLIN + 1);
    long nlin;
    
    check(argc == 3, "usage: dfltovoc <dflin> <outname>");
    dfl = openfile(argv[1], RMODE);
    vfa = openfile(txtcat(argv[2], ".vfa"), WMODE);
    vfs = openfile(txtcat(argv[2], ".vfs"), WMODE);
    vms = openfile(txtcat(argv[2], ".vms"), WMODE); 
    if (vms == NULL) vms = stderr;
    
    define_tabela(tb_pal, "'abcdefghijklmnopqrstuvwxyzáàãâéêíóôõúüç-");
    define_tabela(tb_cat, "ABCDEFGHIJKLMNOPQRSTUVWXYZ?2-.");
    define_tabela(tb_ord, "123456789");
    
    while (le_linha(lin, MAXLLIN, &nlin))
      { if (num_linhas > MAXREG) 
          { fprintf(stderr, "excesso de linhas na entrada\n");
            return(1);
          }
        num_linhas++;
        processa_linha(lin, nlin); 
      }
   
    fprintf(vms, "\n");
    fprintf(vms, "%8ld linhas no dicionario de entrada\n", num_linhas);
    fprintf(vms, "\n");
    fprintf(vms, "%8ld linhas vazias\n", num_vazias);
    fprintf(vms, "%8ld linhas de flexões inexistentes\n", num_inexistentes);
    fprintf(vms, "%8ld linhas malformadas\n", num_erros);
    fprintf(vms, "\n");
    fprintf(vms, "%8ld palavras válidas no vocabulário\n", num_palavras);
    fprintf(vms, "%8ld entradas analítcas gravadas\n", num_analiticas);
    fprintf(vms, "%8ld entradas sintéticas gravadas\n", num_sinteticas);
    return(0);
  } 
  
FILE *openfile(char *name, char *mode)
  { FILE *f; 
    if (vms != NULL) fprintf(vms, "opening %s  mode = %s\n", name, mode);
    if (vms != stderr) fprintf(stderr, "opening %s  mode = %s\n", name, mode);
    f = fopen(name, mode);
    check(f != NULL, txtcat("unable to open ", name));
    return(f);
  }
 
bool_t le_linha(byte *lin, int maxlen, long *nlin)
  { int c = getc(dfl);
    byte *p = lin;
    long n = 0;
    if (c == EOF)
      { check(ferror(dfl) == 0, "erro de leitura em stdin"); 
        *nlin = 0;
        return(false);
      }
    else
      { while((c != EOF) && (c != '\n'))
          { if (n < maxlen) *p = (byte) c; 
            n++; p++;
            c = getc(dfl);
          }
        if (c == EOF)
          { check(ferror(dfl) == 0, "erro de leitura em stdin"); 
            check(false, "ultima linha de stdin não termina com NEWLINE");
          }
        *nlin = n;
        return(true);
      }
  }
  
void mens_linha(char *msg, long num, byte *lin, long nlin)
  { fprintf(vms, "** linha %8ld", num);
    fprintf(vms, ": %s\n", msg);
    if (lin != NULL) imprime_linha(lin, nlin);
    check(num_erros < MAXERR, "excesso de erros");
  }

void processa_linha(byte *lin, long nlin)
  { byte *pal, *ord, *cls, *flx, *can;
    byte *sufpal, *sufcan;
    long npal, nord, ncls, nflx, ncan, nsufpal, nsufcan;
    
    /* Verifica tamanho e caracteres: */
    if (nlin > MAXLLIN) 
      { mens_linha("muito comprida (ignorada)", num_linhas, lin, nlin);
        num_erros++;
        return;
      }
    if (nlin == 0) 
      { mens_linha("vazia (ignorada)", num_linhas, lin, nlin);
        num_vazias++;
        return;
      }

    /* Divide em campos (pal,npal) (cls,ncls),(flx,nflx) (can,ncan) */
    { byte *r = lin; long n = nlin; 
    
      /* Abre parênteses: */
      if (! testa_delim(&r, &n, '(', lin, nlin)) return;
        
      /* Palavra: */
      if (! cata(&r, &n, &pal, &npal, tb_pal, MAXLPAL, lin, nlin)) return; 
      
      /* Separador: */
      if (! testa_delim(&r, &n, ',', lin, nlin)) return;
      
      /* Dígito de ordem */
      if (! cata(&r, &n, &ord, &nord, tb_ord, MAXLORD, lin, nlin)) return; 
      
      /* Separador: */
      if (! testa_delim(&r, &n, ',', lin, nlin)) return;
      
      /* Classe gramatical */
      if(! cata(&r, &n, &cls, &ncls, tb_cat, MAXLCLS, lin, nlin)) return; 
      
      /* Separador: */
      if (! testa_delim(&r, &n, ',', lin, nlin)) return;

      /* Código de flexão */
      if(! cata(&r, &n, &flx, &nflx, tb_cat, MAXLFLX, lin, nlin)) return; 
      
      /* Separador: */
      if (! testa_delim(&r, &n, ',', lin, nlin)) return;

      /* Forma canônica: */
      if(! cata(&r, &n, &can, &ncan, tb_pal, MAXLPAL, lin, nlin)) return; 
      
      /* Fecha parêntese: */
      if (! testa_delim(&r, &n, ')', lin, nlin)) return;

      /* fim de linha: */
      if (n != 0)
        { mens_linha("lixo no fim da linha (ignorada)", num_linhas, lin, nlin);
          num_erros++;
          return;
        }
    }

    /* Campos vazios: */
    if (nord == 0) 
      { mens_linha("ordem vazia (ignorada)", num_linhas, lin, nlin);
        num_erros++;
        return;
      }

    if (ncls == 0) 
      { mens_linha("classe vazia (ignorada)", num_linhas, lin, nlin);
        num_erros++;
        return;
      }

    if (ncan == 0)
      { mens_linha("forma canônica vazia", num_linhas, lin, nlin);
        num_erros++;
        return;
      }
    
    /* Determina maior prefixo comum: */
    sufpal = pal; nsufpal = npal;
    sufcan = can; nsufcan = ncan;
    while((nsufpal > 0) &&  (nsufcan > 0) && ((*sufpal)==(*sufcan)))
      { sufpal++; nsufpal--; sufcan++; nsufcan--; }

    if (nsufpal > 35)
      { mens_linha("sufixo da palavra muito longo", num_linhas, lin, nlin);
        num_erros++;
        return;
      }
    if (nsufcan > 35)
      { mens_linha("sufixo da forma canônica muito longo", num_linhas, lin, nlin);
        num_erros++;
        return;
      }

    if (npal == 0) 
      { mens_linha("entrada de forma inexistente", num_linhas, lin, nlin);
        num_inexistentes++;
      }
    
    /* Grava resultado: */
    num_palavras++;
    fprintf(vfa, 
      "%.*s%s%.*s%.*s%.*s%s%c%.*s#\n", 
      (int)npal, pal, 
      SEPA1, 
      (int)nord, ord, 
      (int)ncls, cls, 
      (int)nflx, flx, 
      SEPA2, 
      digit[nsufpal], 
      (int)nsufcan, sufcan
    );
    num_analiticas++;
    fprintf(vfs, 
      "%.*s%s%.*s%.*s%s%c%.*s#\n", 
      (int)ncan, can, 
      SEPS1, 
      (int)ncls, cls, 
      (int)nflx, flx, 
      SEPS2, 
      digit[nsufcan], 
      (int)nsufpal, sufpal
    );
    num_sinteticas++;
  } 
  
bool_t testa_delim(byte **buf, long *nbuf, char del, byte *lin, long nlin)
  {
    if (((*nbuf) == 0) || (((char)(**buf)) != del))
      { mens_sintaxe("delimitador não encontrado", nlin - (*nbuf), lin, nlin);
        num_erros++;
        return (false);
      } 
    else
      { (*buf)++; (*nbuf)--; return (true); }
  }

bool_t cata(
    byte **buf, long *nbuf, 
    byte **res, long *nres, 
    bool_t tb[], long maxres,
    byte *lin, long nlin
  )
  { (*res) = (*buf);
    (*nres) = 0;
    while (((*nbuf) > 0) && (tb[**buf]))
      { if ((*nres) >= maxres) 
          { mens_sintaxe("campo muito longo (ignorada)", nlin - (*nbuf), lin, nlin);
            num_erros++;
            return (false);
          }
        (*buf)++; (*nbuf)--; (*nres)++; 
      }
    return(true);
  }
    
void mens_sintaxe(char *msg, long col, byte *lin, long nlin)
  { fprintf(vms, "** linha %8ld coluna %3ld", num_linhas, col);
    fprintf(vms, ": %s\n", msg);
    imprime_linha(lin, nlin);
    check(num_erros < MAXERR, "excesso de erros");
  }
  
void imprime_linha(byte *lin, long nlin)
  { long n = nlin; 
    char *ell;
    fprintf(vms, "\"");
    if (n > MAXLLIN) n = MAXLLIN;
    if (n > 40) { n = 40-3; ell = "..."; } else { ell = ""; }
    while (n > 0) 
      { if (((*lin) >= 32) && ((*lin) != 127))
          { fputc((char)(*lin), vms); }
        else
          { fputc((char)(*lin), vms); }
        n--; lin++;
      }
    fprintf(vms, "%s\"\n", ell);
  }
 
void define_tabela(bool_t tb[], char *lets)
  { int k;
    for (k=0; k<256; k++) tb[k] = false;
    while ((*lets)!='\0') {tb[(byte)(*lets)] = true; lets++; }
  }

char *txtcat(char *a, char *b)
  { char *r = malloc(strlen(a) + strlen(b) + 1);
    check(r != NULL, "memory exhausted");
    strcpy(r, a); strcat(r, b); 
    return(r);
  }

  
