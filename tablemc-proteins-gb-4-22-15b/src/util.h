#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED
#include <cstddef>
#include <cstdio>

#define RAD_TO_DEG                  57.2957795//converts radians to degrees
#define DEG_TO_RAD                  1.74532925e-2//converts degrees to radians
#define MB                          1024*1024
#define FALSE 0
#define TRUE 1
#define ALLOCATE(count,type)        ((type *) checkalloc((count),sizeof(type)))

void die(void);
void * checkalloc(size_t count, size_t size);
void * checkrealloc(void * ptr, size_t count, size_t size);
void trim_string(char * s);
void fill_string(char * s,size_t size);
char yesno(int x);
void strlower(char * s);
double atom_mass(char element);
char * read_multiline(FILE * input);
unsigned int digital_crc32(const unsigned char *buf, size_t len);
unsigned long read_random_seed(void);
#if defined(__unix__)
double convtime(struct timeval time);
#endif

#if defined(PARALLEL) || defined(EXCHANGE)
void parallel_init(int * mynod, int * numnod, char * outputfmt);
void parallel_finish();
#endif


#ifdef TIMERS

#define TIMER_NONE                    -1
#define TIMER_MC_MOVE                 0
#define TIMER_NT_BONDS                1
#define TIMER_NT_ANGLES               2
#define TIMER_NT_DIHEDRALS            3
#define TIMER_NT_IMPROPERS            4
#define TIMER_NT_VDW_ELEC             5
#define TIMER_COV_TABLES              6
#define TIMER_PER_ATOM_BORN_RADII     7
#define TIMER_PER_FRAGMENT_BORN_RADII 8
#define TIMER_CONVERT_BORN_RADII      9
#define TIMER_CHECK_CUTOFF            10
#define TIMER_INT_EXACT_PREP          11
#define TIMER_INT_EXACT               12
#define TIMER_INT_PREP                13
#define TIMER_INT_ORIENT              14
#define TIMER_INT_TRANS               15
#define TIMER_INT_INDEX               16
#define TIMER_INT_LOOKUP              17
#define TIMER_INT_OTHER               18
#define TIMER_NB_LIST                 19
#define TIMER_OTHER                   20
#define NTIMERS                       21

struct timer {
    unsigned long long int stop_count;
    unsigned long long int last_started;
    unsigned long long int total_ticks;
};

static struct timer timers[NTIMERS];
static volatile int current_timer;
static const char * timer_names[NTIMERS] = {"MC moves","Bonds","Angles","Dihedrals","Impropers","Non tab. VDW/Elec","Backbone tables",
    "Per-atom Born radii","Per-frag Born radii","Convert Born radii",
    "Check cutoff","Exact prep","Exact interaction","Table prep","Table orientational","Table translational","Table index calc","Table lookup","Other interaction","Nonbond list update","Other"};
static double overhead;

void switch_timer(int timer);
void init_timers(void);
void stop_current_timer(void);
void print_timers(void);
unsigned long long get_rdtsc(void);

#endif //TIMERS
#endif // UTIL_H_INCLUDED
