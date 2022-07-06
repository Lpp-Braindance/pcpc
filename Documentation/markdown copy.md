# N-body

---

Progetto di Programmazione Parallela Concorrente e su Cloud

Nome: Antonello
Cognome: Luppolo

---

## Sommario
- [N-body](#n-body)
  - [Sommario](#sommario)
  - [Introduzione](#introduzione)
  - [Descrizione soluzione proposta](#descrizione-soluzione-proposta)
  - [Dettagli implementativi](#dettagli-implementativi)
    - [Rappresentazione delle particelle](#rappresentazione-delle-particelle)
    - [Distribuzione del carico di lavoro](#distribuzione-del-carico-di-lavoro)
    - [Fase di inizializzazione](#fase-di-inizializzazione)
    - [Fase di simulazione](#fase-di-simulazione)
  - [Istruzioni per l'esecuzione](#istruzioni-per-lesecuzione)
  - [Correttezza del programma](#correttezza-del-programma)
  - [Discussione dei risultati](#discussione-dei-risultati)
    - [Strong scalability](#strong-scalability)
    - [Weak scalability](#weak-scalability)
  - [Conclusioni](#conclusioni)
  
---

## Introduzione

  In generale il problema N-body consiste nel prevedere il moto (posizione e velocità) di un gruppo di N-oggetti che interagiscono indipendentemente l'uno dall'altro, sotto l'influenza di forze fisiche, come ad esempio la  forza di gravità. Molti fenomeni fisici, astrofisici e chimici, possono essere simulati con un sistema di particelle in cui ogni particella interagisce con le altre secondo leggi fisiche.
  Per studiare questi fenomeni vengono utilizzati programmi che simulano il comportamento delle particelle. Tali programmi in:

- **input** prendono un insieme di particelle, dove per ognuna viene specificata la sua posizione nello spazio e la sua velocità all'inizio della simulazione;
- **output** restituiscono la posizione e velocità di ciascuna particella alla fine della simulazione.

Per simulare il comportamento di n corpi sono stati implementati diversi algoritmi con diversi gradi di complessità e approssimazione. Tra questi prendiamo in considerazione il codice realizzato da Mark Harris, reperibile al seguente link: [mini-nbody](https://github.com/harrism/mini-nbody). Si tratta di un programma scritto in linguaggio C che:

- esegue la simulazione su un numero di particelle che può essere specificato dall'utente tramite riga di comando, oppure su un numero predefinito di particelle, nel caso in cui non venisse fornito alcun input.
- tramite un algoritmo pseudocasuale vengono inizializzati i campi di ogni particella
- la simulazione prevede 10 iterazioni dove ad ognuna:
  - si tiene traccia del tempo di inizio dell'iterazione
  - in maniera sequenziale ogni particella viene fatta interagire con tutte le altre per calcolare la forza e la posizione che assumerà all'iterazione successiva
  - viene stampato sullo standard output il tempo impiegato per il completamento dell'iterazione i-esima
- in fine viene stampato, sempre su standard output, il rapporto iterazioni/secondo.

La soluzione fornita da tale programma è quadratica nel numero di particelle.
È possibile ottimizzare tale codice utilizzando la libreria OpenMPI, in modo da parallelizzare il calcolo.

## Descrizione soluzione proposta

Di seguito viene mostrata una versione ottimizzata del programma sopra citato che fa uso della libreria OpenMPI, mediante la quale si va a distribuire il calcolo della forza delle particelle tra n processi nel seguente modo:

- Innanzi tutto i processi comunicano tra loro attraverso il communicator MPI_COMM_WORLD. All'interno del gruppo associato a questo communicator i processi vengono distinti in base al loro rank in due "categorie":
  - MASTER: un solo processo appartiene a questa categoria e viene identificato con rank = 0;
  - SLAVE: tutti gli altri invece appartengono a quest'altra categoria ed hanno rank > 0;
- tutti i processi contribuiscono al calcolo della forza delle particelle.
- la simulazione viene fatta su un numero preciso di particelle se tale valore viene specificato dall'utente tramite riga di comando. Altrimenti viene fatta su un numero predefinito di particelle, indicato all'interno del programma.
- ciascun processo, una volta a conoscenza del numero di particelle e processi, esegue il calcolo per la distribuzione del carico di lavoro. Cioè calcola la propria porzione di particelle e quella degli altri processi.
In seguito procede con la simulazione, la quale prevede 10 iterazioni, dove in ciascuna ogni processo:  
  1) tramite invocazioni della funzione collettiva non bloccante MPI_Ibcast invia solo le "coordinate"(valori : x, y, z) delle sue particelle a tutti i processi e riceve le altre da essi. Vengono scambiate solo le coordinate delle particelle poiché ciascun processo per il calcolo della velocità (valori vx, vy, vz) e della forza di ciacuna sua particella ha bisogno solo delle coordinate di tutte le altre.
  2) durante la fase di comunicazione inizia il calcolo della forza relativo alla sua porzione di particelle. Ovvero fa interagire ciascuna sua particella con tutte le altre al suo interno.
  3) Resta in attesa di ricevere una o più porzioni tramite la funzione MPI_Waitany. Non appena riceve una o più porzioni riprende il calcolo della forza, utilizzando le particelle appena ricevute, facendole interagire con le proprie. Ripete il passo 3 fin quando ci sono ancora richieste di ricezione da completare.
  4) terminate le fasi di ricezione e calcolo della forza, il processo può proseguire con l'aggiornamento dei valori(coordinate e velocità) della propria porzione di particelle. Tali risultati sono necessari per il prossimo passo della simulazione. Quindi come prima cosa si vanno ad aggiornare i valori relativi alla velocità e i processi SLAVE inviano tali risultati al MASTER tramite la funzione non bloccante MPI_Igatherv. Solo il MASTER ha bisogno di tutti i valori delle particelle poichè esso deve collezionarili per poi stamparli su un file. C'è da tener presente che gli SLAVE dalla seconda iterazione in poi, per aggiornare i propri valori devono attendere (tramite la MPI_Wait) che la richiesta di invio precedente(relativa alla MPI_Igatherv) venga completata.
  
## Dettagli implementativi

### Rappresentazione delle particelle

Le particelle nel programma sequenziale vengono rappresentate tramite la struttura :

```C
typedef struct { float x, y, z, vx, vy, vz; } Body;
```

Mentre nella versione parallela del programma, vengono rappresentate attraverso due *struct*:

```C

typedef struct { float x, y, z; } BodyPosition;
typedef struct { float vx, vy, vz; } BodyVelocity;

```

La prima contiene le informazioni riguardanti le coordinate nello spazio e l'altra invece i valori relativi alla velocità. Dunque per le particelle vengono allocati due array, uno per le coordinate e l'altro per le velocità(entrambi con dimensione pari al numero di particelle dato in input). Questa suddivisione è stata fatta perchè il processo per calcolare la forza e la velocità delle proprie particelle ha bisogno di conoscere solo le coordinate di tutte le altre particelle. L'unico processo che invece ha bisogno di conoscere otlre alle coordinate anche gli altri valori è il MASTER, poiché esso deve collezionare tutti i risultati per poi stamparli su un file. Grazie a queste osservazioni è stato possibile ridurre considerevolmente l'overhead di comunicazione allo stretto necessario.
Visto che ogni particella nella sua interezza è una sestupla di float, dove ciascun campo rappresenta. Dunque nel caso del programma sequenziale, per quanto detto è possibile dichiarare l'insieme delle particelle nel seguente modo: 

```C
  Body p[nBodies];
```

mentre nella versione parallela del programma oltre a dichiarare i due array per le coordiante e le velocità, vengono dichiarati i corrispondenti tipi di dati derivati MPI per utilizzare queste strutture nelle funzioni di comunicazione MPI:

```C
  MPI_Datatype BODY_POS, BODY_VEL;
  MPI_Datatype strcut_types[1];
  int blockcounts[1];
  MPI_Aint offsets[1], lb, extent;
  offsets[0] = 0;
  strcut_types[0] = MPI_FLOAT;
  blockcounts[0] = 3;
  
  MPI_Type_create_struct(1, blockcounts, offsets, strcut_types, &BODY_POS);
  MPI_Type_commit(&BODY_POS);

  MPI_Type_create_struct(1, blockcounts, offsets, strcut_types, &BODY_VEL);
  MPI_Type_commit(&BODY_VEL);

  BodyPosition body_pos[nBodies];
  BodyVelocity body_vel[nBodies];
  
```

### Distribuzione del carico di lavoro

```C
  ...
  // CALCULATE WORKLOAD DISTRIBUCTION
  int portion = nBodies / n_workers;
  int rest = nBodies % n_workers;
  int portions_sizes[n_workers], portions_starts[n_workers];
  calculatePortions(portions_sizes, portions_starts, n_workers, portion, rest);
  ...
```

Indichiamo con *n_workers* il numero di processi che contribuiscono al calcolo della forza delle particelle. 

Come già accennato, ogni processo contribuisce alla computazione e il processo MASTER, in più rispetto agli altri tiene traccia del tempo impiegato dall'iterazione i-esima della simulazione e colleziona i risultati per poi stamparli su un file.

 a video tali informazioni. Il motivo per il quale è stata fatta tale scelte deriva dal fatto che: il tempo necessario per l'iterazione i-esima, in cui il MASTER oltre a partecipare al calcolo della forza delle particelle esegue anche le operazioni appena dette, è minore in confronto ad una nella quale il MASTER non contribuisce al calcolo della forza, ma riceve solo i risultati e tiene solo traccia del tempo facendone la stampa a video. Poichè gli altri n-1 processi avrebbero più carico di lavoro e il processo MASTER resterebbe in stato di hidle fino al ricevimento dei risultati. Dunque il tempo complessivo del programma in cui il MASTER "non lavora" è maggiore rispetto a quando esso contribuisce al calcolo della forza. Ecco perchè in questo caso specifico conviene che tutti i processi contribuiscano al calcolo.


Il numero di particelle viene diviso tra gli n processi andando a dividere la taglia dell'input per il numero di n_workers. Calcolata la porzione che dovrà essere assegnata a ciascun processo e l'eventuale resto, vengono allocati gli array *procs_portions_sizes[ ]* e *procs_portions_starts[ ]*. Questi array servono per tener traccia delle porzioni di ciascun processo per la fase di invio e ricezione delle altre particelle(questo aspetto verrà approfondito in più avanti, quando verrà trattata la fase di comunicazione tramite chiamate collettive). Dopodichè con l'invocazione della funzione *calculatePortion()* viene fatto il calcolo delle porzioni da assegnare a ciascun processo come mostrato di seguito

```C
void calculatePortions(int procs_portions_sizes[], int procs_portions_starts[], int n_workers, int portion, int rest)
{
    for (int i = 0; i < n_workers; i++)
    {
        procs_portions_sizes[i] = portion;
        if (rest > 0)
        {
            procs_portions_sizes[i]++;
            rest--;
        }
    }

    procs_portions_starts[0] = 0;
    for (int i = 1; i < n_workers; i++)
        procs_portions_starts[i] = procs_portions_starts[i - 1] + procs_portions_sizes[i - 1];
}
```

Se il numero di particelle non è divisibile per il numero di processi allora il resto delle particelle viene distribuito a partire dal 1° processo (il cui rank = 0) in poi. Cioè viene assegnata una particella in più, delle restanti, a tali processi.  Altrimenti le particelle vengono distribuite equamente tra gli n processi.
Dunque per ogni processo viene calcolata la propria porzione e l'indice dal quale inizia il suo intervallo nel corrisponte array di particelle (descritto al punto precedente).

### Fase di inizializzazione

```C
    typedef struct 
    { 
        int rank, own_portion, start_own_portion;
        int *procs_portions_sizes, *procs_portions_starts;
        float *Fx, *Fy, *Fz;

    } ProcVariables;

    ...

    // DEFINE PROCES VARIABLES
    ProcVariables proc;
    proc.rank = worker_rank;
    proc.own_portion = portions_sizes[worker_rank];
    proc.start_own_portion = portions_starts[worker_rank];
    float Fx[proc.own_portion], Fy[proc.own_portion], Fz[proc.own_portion];
    proc.Fx = Fx;  proc.Fy = Fy;  proc.Fz = Fz;
    proc.procs_portions_sizes = portions_sizes;
    proc.procs_portions_starts = portions_starts;
    // INIT OWN BODY PORTION
    determisticInitBodiesSplit(body_pos, body_vel, proc.own_portion, proc.start_own_portion);
    
    ...


```

Viene definita la variabile proc, la quale è una struttura usata per raggruppare le variabili che servono al processo in esecuzione come ad esempio la taglia della sua porzione, l'indice di inizio del suo intervallo di particelle, e gli array *Fx*, *Fy*, *Fz* che servono per tener traccia dei valori intermedi relativi al calcolo della forza delle proprie particelle, durante i vari calcoli(tra le diverse invocazioni della funzione *bodyForceSplit*) con le altre particelle che vengono ricevute dagli altri processi(l'utilità di tali variabili verrà approfondita più avanti).
Ogni processo, una volta definito il proprio intervallo, inizializza le proprie particelle attraverso un algoritmo deterministico, ovvero tramite la funzione *determisticInitBodiesSplit*. A differenza del programma di partenza, che fa uso di un algoritmo che assegna valori alle particelle in maniera casuale, è stato utilizzato un algoritmo deterministico per poter valutare la correttezza del programma. Nel dettaglio la funzione utilizzata a tale scopo è la seguente:

```C

void determisticInitBodiesSplit(BodyPosition *body_pos, BodyVelocity *body_vel, int own_portion, int start_own_portion)
{
    for (int i = 0; i < own_portion; i++)
    {
        body_pos[start_own_portion + i].x = start_own_portion + i + 1;
        body_pos[start_own_portion + i].y = start_own_portion + i + 1;
        body_pos[start_own_portion + i].z = start_own_portion + i + 1;

        body_vel[start_own_portion + i].vx = start_own_portion + i + 1;
        body_vel[start_own_portion + i].vy = start_own_portion + i + 1;
        body_vel[start_own_portion + i].vz = start_own_portion + i + 1;
    }
}
```

### Fase di simulazione

```C
  // DEFINE REQUESTS FOR PORTIONS EXCHANGES
  MPI_Request bcast_reqs[n_workers];
  MPI_Request body_vel_gath = MPI_REQUEST_NULL;
  int reqs_ranks[n_workers];
  // START i-th PROCES WORK
  for (int iter = 1; iter <= nIters + 1; iter++)
  {
      for (int i = 0; i < n_workers; i++)
          MPI_Ibcast(&body_pos[portions_starts[i]], portions_sizes[i], BODY_POS, i, MPI_COMM_WORLD, &bcast_reqs[i]);

      if (iter > 1)
      {
          for (int i = 0; i < proc.own_portion; i++) {  Fx[i] = 0.0f;  Fy[i] = 0.0f;  Fz[i] = 0.0f;  }
          waitSomeWork2(bcast_reqs, reqs_ranks, n_workers, body_pos, proc);

          if (worker_rank == MASTER && iter > 2 && print_res == 1)
              printResults2(body_pos, body_vel, nBodies);

          integrateVelocitySplit2( body_vel, dt, proc.own_portion, proc.start_own_portion, Fx, Fy, Fz);
          integratePositionSplit2( body_pos, body_vel, dt, proc.own_portion, proc.start_own_portion);
      }
      MPI_Wait(&body_vel_gath, MPI_STATUS_IGNORE);
      MPI_Igatherv(&body_vel[proc.start_own_portion], proc.own_portion, BODY_VEL, &body_vel[proc.start_own_portion], portions_sizes, portions_starts, BODY_VEL, MASTER, MPI_COMM_WORLD, &body_vel_gath); 
      MPI_Barrier(MPI_COMM_WORLD);

  }
  // PROCESES SEND ONLY LAST RESULT TO THE MASTER
  MPI_Gatherv(&body_pos[proc.start_own_portion], proc.own_portion, BODY_POS, &body_pos[proc.start_own_portion], portions_sizes, portions_starts, BODY_POS, MASTER, MPI_COMM_WORLD);
  MPI_Gatherv(&body_vel[proc.start_own_portion], proc.own_portion, BODY_VEL, &body_vel[proc.start_own_portion], portions_sizes, portions_starts, BODY_VEL, MASTER, MPI_COMM_WORLD);

  if (worker_rank == MASTER && print_res == 1)
      printResults2(body_pos, body_vel, nBodies); // stampiamo qui fuori perchè quando usciamo dal ciclo for la igather sta ancora lavorando e quindi qui dobbiamo stampare i risultati
}

```

In questa fase ongi processo invia le coordinate della propria porzione di particelle a tutti gli altri, utilizzando la funzione collettiva non bloccante MPI_Ibcast e tiene traccia delle corrispondenti richieste. Nel mentre inzia a calcolare la forza delle sue particelle. Ovvero come prima cosa va ad inizializzare gli array Fx, Fy, Fz, il cui scopo è di memorizzare i "valori intemedi" della forza delle proprie particelle durante i vari calcoli (tra le varie invocazioni della funzione *bodyForceSplit* come verrà a breve) con le altre che vengono ricevute da ciascun processo. Dunque Fx, Fy, Fz si prebbero pensare come una sorta di "accumulatori temporanei" corrispondenti a ognuna particella per calcolare la sua forza ad ogni iterazione della simulazione. In seguito il processo è pronto ad iniziare il calcolo vero e proprio della forza di ciascuna sua particella invocando la funzione *waitSomeWork*:




tramite invocazioni della funzione collettiva non bloccante MPI_Ibcast invia solo le "coordinate"(valori : x, y, z) delle sue particelle a tutti i processi e riceve le altre da essi.

Qui ogni processo tramite invocazioni della funzione collettiva non bloccante MPI_Ibcast invia la propria porzione di particelle a tutti gli altri processi, e riceve da essi le altre porzioni. In questa fase si tiene traccia delle corrispondenti richieste.
Durante la fase di comunicazione il processo inizializza gli array Fx, Fy, Fz settando tutti i valori a 0. Lo scopo di tali array è memorizzare i "valori intermedi" della forza delle proprie particelle durante i vari calcoli (tra le varie invocazioni della funzione bodyForceSplit) con le altre che vengono ricevute. Dunque Fx, Fy, Fz si potrebbero pensare come a una sorta di "accumulatori temporanei" corrispondenti a ciascuna particella della porzione del processo, per calcolare la loro forza ad ogni iterazione. Fatto ciò, il processo è pronto ad iniziare il calcolo della forza della propria porzione invocando la funzione bodyForceSplit:

```C
void bodyForceSplit(Body *p, float dt, int own_portion, int start_own_portion, int start, int end, float Fx[], float Fy[], float Fz[])
{
  for (int i = 0; i < own_portion; i++)
  {
    for (int j = start; j < end; j++) // start = inizio nuova porzione da confrontare end = fine porzione da confrontare
    {
      float dx = p[j].x - p[start_own_portion + i].x;
      float dy = p[j].y - p[start_own_portion + i].y;
      float dz = p[j].z - p[start_own_portion + i].z;

      float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
      float invDist = 1.0f / sqrt(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx[i] += dx * invDist3;
      Fy[i] += dy * invDist3;
      Fz[i] += dz * invDist3;
    }
  }
}
```

La funzione bodyForeceSplit prende in input essenzialmente l'array di particelle tramite il parametro p (puntatore che punta all'array *buf*), un valore costante dt, la taglia della porzione del processo i-esimo (numero di particelle della sua porzione), l'indice dal quale inizia la sua porzione nell'array di particelle, l'inizio(parametro *start*) e la fine(parametro *end*) dell'intervallo (in *buf*) della nuova porzione da utilizzare per il calcolo, ed infine gli array Fx, Fy, Fz che contengono i valori intermedi.
Inizialmente il processo i-esimo invoca la funzione bodyForceSplit passando come parametri *start* ed *end* il proprio intervallo, ovvero i suoi indici di inizio e fine in *buf*. Dunque comincia a calcolare la forza di ciascuna particella della sua porzione, facendo "interagire" ognuna di esse con tutte quante le altre al suo interno. Dopodiché tramite la funzione wait_some_work , sotto riportata:

```C
void waitSomeWork(MPI_Request bcast_recv[], int request_rank_indices[], int requests_ranks[], int n_recv_req, int rank, Body *p, float Fx[], float Fy[], float Fz[], int proc_portion_size[], int proc_portion_start[])
{
  int ready_req = 0;
  int count = n_recv_req;
  int start, end;
  int req_compl;
  const float dt = 0.01f;
  int own_portion = proc_portion_size[rank];
  int start_own_portion = proc_portion_start[rank];

  MPI_Waitsome(n_recv_req, bcast_recv, &ready_req, request_rank_indices, MPI_STATUS_IGNORE);
  while (ready_req != MPI_UNDEFINED)
  {
      count -= ready_req;
      for (int i = 0; i < ready_req; i++)
      {
          req_compl = requests_ranks[request_rank_indices[i]];
          start = proc_portion_start[req_compl];
          end = start + proc_portion_size[req_compl];

          bodyForceSplit(p, dt, own_portion, start_own_portion, start, end, Fx, Fy, Fz);
      }
      MPI_Waitsome(n_recv_req, bcast_recv, &ready_req, request_rank_indices, MPI_STATUS_IGNORE);
  }
  
  for (int i = 0; i < count; i++)
  {
      req_compl = requests_ranks[request_rank_indices[i]];
      start = proc_portion_start[req_compl];
      end = start + proc_portion_size[req_compl];

      bodyForceSplit(p, dt, own_portion, start_own_portion, start, end, Fx, Fy, Fz);
  }

}
```

attendere il completamento delle richieste di ricezione delle altre porzioni degli n-1 processi per completare il calcolo della forza. Il primo motivo per il quale si tiene traccia delle richieste di comunicazione è per far in modo che il calcolo continui man mano che le particelle vengono ricevute. A tal proposito viene utilizzata la funzione MPI_Waitsome, la quale attende che almeno una delle operazioni associate agli handle attivi, nell'array di richieste(*bcast_recv*), sia stata completata e restituisce nella variabile *ready_req* il numero di richieste completate. Inoltre va a scrivere nell'array di indici (*request_rank_indices*), nelle posizioni da 0 a *ready_req*, i corrispondenti indici delle richieste completate. Per risalire al rank del worker relativo a ciascuna richiesta, ogni processo, prima di iniziare la simulazione, inizializza l'array requests_ranks.
Quando una richiesta è stata completata viene deallocata, e l'handle associato viene impostato a MPI_REQUEST_NULL. Se l'array di richieste (*bcast_recv*) non contiene handle attivi, la chiamata ritorna immediatamente restituendo in *ready_req* il valore MPI_UNDEFINED. Dunque la prima volta che viene chiamata MPI_Waitsome se ci sono ancora richieste che non sono state ancora completate si entra nel ciclo while. All'interno del ciclo viene decrementata la variabile *count*, la quale tiene traccia delle richieste che devono essere ancora completate. *count* ha un ruolo importante che a breve verrà chiarito. Poi si va a scorrere l'array degli indici delle richieste completate, dove per ognuna di queste, si va a calcolare l'intervallo di particelle (in *buf*) del processo corrispondente, sul quale andare ad invocare la funzione bodyForceSplit per riprendere il calcolo della forza delle particelle del processo in esecuzione. Terminato il calcolo intermedio utilizzando le particelle ricevute, viene invocata nuovamente la funzione MPI_Waitsome e tutto ciò viene iterato fin quando non ci sono più richieste attive e quindi fin quando l'ultima MPI_Waitsome restituisce in *ready_req* il valore MPI_UNDEFINED. Ma quando si esce dal *while* non sono state ancora utilizzate le ultime particelle ricevute! È qui che si nota l'utilità della variabile *count*, la quale ha un duplice scopo. Il primo è quello di andare appunto a completare il calcolo della forza utilizzando le ultime particelle ricevute tramite il ciclo *for*, il quale non fa altro che scorrere l'array di richieste completate e invocare la funzione bodyForceSplit sull'intervallo del processo corrispondente. L'altro motivo è che se la prima volta che viene invocata la MPI_Waitsome (prima di entrare nel *while*), se tutte le richieste sono state completate, il valore restituito in *ready_req* sarà MPI_UNDEFINED. Per cui non si entrerà nel ciclo *while* e di conseguenza le particelle ricevute "associate a queste richieste" non verrebbero utilizzate. Ecco perché all'inizio il valore *count* viene posto uguale al numero di richieste da aspettare, il quale viene passato in input alla funzione(parametro *n_recv_req*). Così facendo tramite lo stesso ragionamento di prima e cioè tramite lo stesso ciclo *for*, sopra descritto, vengono utilizzate le particelle ricevute.
Una volta che il processo i-esimo ha ricevuto tutte le restanti porzioni ed ha terminato il calcolo della forza delle proprie particelle, esce dalla funzione waitSomeWork ed è pronto per aggiornare i loro valori(posizioni e forza) tramite la funzione integratePositionSplit: 

```C
void integratePositionSplit(Body *p, float dt, int own_portion, int start_own_portion, float Fx[], float Fy[], float Fz[])
{
    for (int i = 0; i < own_portion; i++)
    {
        p[start_own_portion + i].vx += dt * Fx[i];
        p[start_own_portion + i].vy += dt * Fy[i];
        p[start_own_portion + i].vz += dt * Fz[i];

        p[start_own_portion + i].x += p[start_own_portion + i].vx * dt;
        p[start_own_portion + i].y += p[start_own_portion + i].vy * dt;
        p[start_own_portion + i].z += p[start_own_portion + i].vz * dt;
    }
}
```

A questo punto può essere ripreso il discorso sull'importanza del tener traccia delle richieste di comunicazione. Come preannunciato, la richiesta di invio del processo i-esimo ha un ruolo fondamentale. Il motivo per il quale tale richiesta ha particolare importanza è per far in modo che il processo dopo aver completato il calcolo della forza delle sue particelle, nel caso in cui stesse ancora trasmettendo la sua porzione, non debba aspettare il suo completamento per poi passare all'iterazione successiva. Siccome, arrivati a questo punto, il processo disporrebbe già della porzione necessaria sulla quale lavorare. Per cui, sia il MASTER sia gli SLAVE, prima di aggiornare i loro valori controllano tramite la funzione MPI_Test lo stato della richiesta. Se la richiesta di invio "corrente" non è stata completata si scrivono i risultati in un secondo buffer, il quale verrà poi utilizzato nella prossima iterazione. Così facendo il processo può proseguire con il lavoro evitando inutili attese. Per implementare tale meccanismo, quando si passa all’iterazione successiva si invia la nuova porzione utilizzando un'altra *MPI_Request*(diversa dalla precedente). Dunque vengono dichiarate due *MPI_Request* (*bcast_send_prec*, *bcast_send_next*) e un due puntatori a *MPI_Request* (*bcast_pointer_next_req* e *bcast_pointer_prec_req*). Il puntatore *bcast_pointer_next_req* viene utilizzato nella funzione *MPI_Ibcast* corrispondente alla chiamata in cui il processo in esecuzione invia la sua porzione. C’è da tener presente che dalla 2° iterazione in poi, se sia la richiesta di invio "corrente" che quella "precedente" non sono state ancora completate, in questo caso il processo i-esimo dovrà aspettare il completamento della richiesta precedente, visto che non può scrivere nell’altro buffer.
Al termine di ogni iterazione in cui si vanno a scrivere i risultati nell'altro buffer, per garantire la relazione di precedenza sia per le richieste sia per i buffer di invio/ricezione viene fatto lo scambio dei rispettivi puntatori.
Altrimenti, in tutti gli altri casi, se il processo una volta completato il calcolo della forza delle sua porzione, ha già terminato la fase di trasmissione, viene riutilizzato il buffer corrente anche per l'iterazione successiva.

## Istruzioni per l'esecuzione

compilazione :

```bash
mpicc nbody.c -lm -o {nome eseguibile}
es: mpicc nbody.c -lm -o nbody_split.out
```

esecuzione senza output su file:

```bash
mpirun -np {numero processi} {nome eseguibile} {numero particelle}
es: mpirun -np 4 nbody_split.out 40000
```

esecuzione con output scritto su file per testare la correttezza del programma:

```bash
mpirun -np {numero processi} {nome eseguibile} {numero particelle} -t
es: mpirun -np 4 nbody_split.out 40000 -t
```

## Correttezza del programma

La correttezza del programma viene verificata attraverso n esecuzioni dello stesso. Ogni esecuzione viene fatta utilizzando il comando

```bash
mpirun -np {numero processi} {nome eseguibile} {numero particelle} -t
```

Dove ad ogni esecuzione viene dato in input al programma la stessa quantità di particelle e il corrispondente numero di p processi da utilizzare(assegnato in maniera progressiva da 1 a n). Quando viene specificato nel comando il parametro -t, all'inizio dell'esecuzione il programma genera un file nel quale poi vengono scritti i risultati dalla simulazione. Il nome con il quale viene creato il file dipende dal numero di processi con cui si esegue il programma. Per semplicità i file vengono nominati parallel_{numero processi} es:

```bash
es: mpirun -np 4 nbody_split.out 40000 -t
genera il file parallel_4.txt
```

Per verificare che l'esecuzione fatta con uno o più processi generi lo stesso output, al termine di tutte le esecuzioni, il file generato dal programma sequenziale (parallel_1.) viene confrontato con tutti gli altri.
Per automatizzare la verifica della correttezza è stato realizzato il seguente script bash:

```bash
#!/bin/bash
nBodies=$1 #prende in input come primo parametro il numero di bodies 
n_proc=$2  #prende in input come secondo parametro il numero di massimo di n processi
#controlla se viene specificato il parametro -t per eseguire il programma che salva i risultati su file
if [[ "$#" -eq 3 && $3 == "-t" ]] 
then
    mpirun -np 1 nbody_split.out $nBodies -t
fi
#per non mostare l'output del programma, ma solo l'output dei test, specificare -dn
if [[ "$#" -eq 4 && $4 == "-dn" ]] 
then
    mpirun -np 1 nbody_split.out $nBodies -t 1>/dev/null
fi
# attendiamo che il programma sequenziale termini la sua esecuzione
wait
printf "sequential program ready\n\n"
# settiamo la variabile proc la quale va a specificare il numero di processi per l'i-esima esecuzione 
proc=2
# testiamo la correttezza del programma eseguendo lo stesso con un numero progressivo di processi 
# cioè eseguiamo il programma con 1,2,3,...,max_nproc 
while [[ ((proc -le n_proc)) ]]
do  
 if [ "$#" -eq 4 ] 
  then
    mpirun -np $proc nbody_split.out $nBodies -t 1>/dev/null
  else
    mpirun -np $proc nbody_split.out $nBodies -t
  fi
  # attendiamo che il programma termini la sua esecuzione 
  wait
  # confrontiamo il risultati che ha generato e memorizzato nel suo corrispondente file con quelli
  # generati dal file sequenziale tramite il comando diff nel quale specifichiamo il parametro -q
  # che indica a diff di fornire l'output solo se ci sono differenze tra i due file
  DIFF=$(diff -q parallel_$proc sequential_program )
  if [ "$DIFF" != "" ] 
  then
    printf "parallel program np %d \!\!\! ERRORE \!\!\! \n\n" $proc
  else
    printf "parallel program np %d --> OK\n\n" $proc
  fi
  ((proc++))
done  

```

## Discussione dei risultati

Le prestazioni della soluzione proposta sono state valute tenendo in considerazione sia la scalabilità forte che quella debole. I benchmark sono stati svolti su un cluster Google Cloud con 6 istanze di macchine virtuali (e2-standard-4), aventi ciascuna 2 core(4vCPU).
Sia per la scalabilità forte che per quella debole sono stati eseguiti 12 esperimenti.

### Strong scalability

La strong scalability viene valutata, mantenendo la taglia dell'input del problema fissa e facendo aumentare progressivamente il numero di processi, con i quali eseguire il programma. Gli esperimenti sono stati condotti su un input fissato di 30000 particelle. Di seguito vengono mostrati i risultati ottenuti:

<div align="center">
<table>
  <tr>
    <td><p weight="400px">

vCPU   |  Input  | Tempo
------ | ------  | -----
  1    |  30000  | 173,624
  2    |  30000  | 123,140
  3    |  30000  | 84,922
  4    |  30000  | 61,565
  5    |  30000  | 53,424
  6    |  30000  | 48,457
  7    |  30000  | 43,463
  8    |  30000  | 37,239
  9    |  30000  | 35,890
  10   |  30000  | 30,699
  11   |  30000  | 22,441
  12   |  30000  | 20,770
  </p></td>
  <td>

  ![alt](./Documentation/strong_scalability_graph.png)

  </td>
    </tr>
</table>
</div>

Dagli esiti si può riscontrare che all'aumentare del numero di processi, il tempo di esecuzione diminuisce fino ad un certo punto. Da lì in poi la velocità di esecuzione del programma aumenta sempre meno. Questo perchè con l'aumentare del numero di processi, il lavoro svolto da ciascuno di essi diminuisce, ma di contro aumenta l'overhead di comunicazione. Dunque con l'aumentare del numero di processi ognuno di essi spende sempre meno tempo per la computazione e sempre più tempo per la comunicazione che è la parte più onerosa dell'esecuzione e che quindi richiede più tempo.

### Weak scalability

La weak scalability invece, viene valutata facendo variare la taglia dell'input al variare del numero di processi facendo in modo che ogni processo abbia sempre lo stesso workload(ovvero la stessa quantità di lavoro). Gli esperimenti sono stati fatti in modo che ogni processo lavorasse su 2500 particelle. Dunque partendo con un processo con 2500 particelle in input, fino ad arrivare a 12 processi con 30000 particelle come input. Di seguito vengono mostrati i risultati ottenuti:

<div align="center">
<table>
  <tr>
    <td><p weight="400px">

vCPU   |  Input  | Tempo
------ | ------  | -----
  1    |  2500   | 1,203
  2    |  5000   | 3,434
  3    |  7500   | 5,130
  4    |  10000  | 6,888
  5    |  12500  | 8,577
  6    |  15000  | 10,413
  7    |  17500  | 12,020
  8    |  20000  | 13,807
  9    |  22500  | 15,519
  10   |  25000  | 17,228
  11   |  27500  | 18,849
  12   |  30000  | 20,592
  </p></td>
  <td>

  ![alt](./Documentation/weak_scalability_graph.png)

  </td>
    </tr>
</table>
</div>

Possiamo notare dai risultati della weak scalability che le prestazioni degradando costantemente di quasi la metà ogni volta che taglia e numero di processi crescono progressivamente. Il problema in questo caso è che man mano che si scala l'overed di comunicazione per ciascun processo aumenta. Più nel dettaglio, il lavoro di computazione per ogni processo rimane lo stesso, ma a questo va aggiunto il tempo di comunicazione, il quale aumenta man mano che si scala, poichè ogni processo deve comunicare con più processi.

## Conclusioni

Come si è potuto riscontrare dalla strong scalability i risultati ottenuti sono soddisfacenti. In termini di tempo di esecuzione sono stati riscontrati miglioramenti. Mentre gli esiti dalla weak scalability non sono stati così tanto positivi visto che vi è un costante incremento, sempre in termini di tempo di esecuzione, di quasi il doppio man mano che si scala.
In conclusione i risultati ottenuti si possono ritenere accettabili vista la complessità dell'algoritmo e in particolare il suo considerevole overhead di comunicazione, il quale è causato dal fatto che ogni processo ha bisogno delle porzioni di tutti quanti gli altri per svolgere il proprio lavoro.
