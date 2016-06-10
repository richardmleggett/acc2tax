/*----------------------------------------------------------------------*
 * File:    acc2tax.c
 * Author:  Richard Leggett (richard.leggett@tgac.ac.uk)
 * Purpose: Convert GI or accession to taxonomy
 * Created: 10 Jan 2013
 * Update:  10 Jun 2016
 *----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h> 

/*----------------------------------------------------------------------*
 * Constants
 *----------------------------------------------------------------------*/
#define MAX_GI 1050000000
#define MAX_NAMES 2000000
#define MAX_PATH 10000
#define VERSION "v0.4"

/*----------------------------------------------------------------------*
 * Globals
 *----------------------------------------------------------------------*/
char database_dir[MAX_PATH];
char input_filename[MAX_PATH];
char output_filename[MAX_PATH];
int is_nucleotide = 1;
int is_protein = 0;
int is_accession = 1;
int is_gi = 0;
long int memory_required = 0;
unsigned int* gi_to_node;
char** names;
unsigned int* nodes;
unsigned int max_gi = MAX_GI;
FILE* acc_fp;
long int acc_file_size;

/*----------------------------------------------------------------------*
 * Function:   usage
 * Purpose:    Report program usage.
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void usage(void)
{
    printf("Bugs/comments: richard.leggett@tgac.ac.uk\n" \
           "\nProvide batch taxonomy information for Genbank IDs or Accessions.\n" \
           "\nOptions:\n" \
           "    [-h | --help]       This help screen.\n" \
           "    [-a | --accession]  Query is accession IDs [default].\n" \
           "    [-d | --database]   Directory containing NCBI taxonomy files.\n" \
           "    [-e | --entries]    Max GI entries (default 1050000000).\n" \
           "    [-g | --gi]         Query is Genbank IDs.\n" \
           "    [-i | --input]      File of IDs (GI or Accession), one per line.\n" \
           "    [-n | --nucleotide] Query IDs are nucleotide [default].\n" \
           "    [-o | --output]     Filename of output file.\n" \
           "    [-p | --protein]    Query IDs are protein.\n" \
           "\n");
}

/*----------------------------------------------------------------------*
 * Function:   parse_command_line
 * Purpose:    Parse command line options
 * Parameters: argc = number of arguments
 *             argv -> array of arguments
 * Returns:    None
 *----------------------------------------------------------------------*/
void parse_command_line(int argc, char* argv[])
{
    static struct option long_options[] = {
        {"accession", no_argument, NULL, 'a'},
        {"database", required_argument, NULL, 'd'},
        {"entries", required_argument, NULL, 'e'},
        {"gi", no_argument, NULL, 'g'},
        {"help", no_argument, NULL, 'h'},
        {"input", required_argument, NULL, 'i'},
        {"nucleotide", no_argument, NULL, 'n'},
        {"output", required_argument, NULL, 'o'},
        {"protein", no_argument, NULL, 'p'},
        {0, 0, 0, 0}
    };
    int opt;
    int longopt_index;
    
    input_filename[0] = 0;
    output_filename[0] = 0;
    database_dir[0] = 0;
    
    while ((opt = getopt_long(argc, argv, "ad:e:ghi:no:p", long_options, &longopt_index)) > 0)
    {
        switch(opt) {
            case 'h':
                usage();
                exit(0);
                break;
            case 'a':
                is_accession = 1;
                is_gi = 0;
                break;
            case 'g':
                is_accession = 0;
                is_gi = 1;
                break;
            case 'i':
                if (optarg==NULL) {
                    printf("Error: Option requires an argument.\n");
                    exit(1);
                }
                strcpy(input_filename, optarg);
                break;
            case 'o':
                if (optarg==NULL) {
                    printf("Error: Option requires an argument.\n");
                    exit(1);
                }
                strcpy(output_filename, optarg);
                break;
            case 'd':
                if (optarg==NULL) {
                    printf("Error: Option requires an argument.\n");
                    exit(1);
                }
                strcpy(database_dir, optarg);
                break;
            case 'e':
                if (optarg==NULL) {
                    printf("Error: Option requires an argument.\n");
                    exit(1);
                }
                max_gi=atoi(optarg);
                break;
            case 'n':
                is_nucleotide = 1;
                is_protein = 0;
                break;
            case 'p':
                is_nucleotide = 0;
                is_protein = 1;
                break;
        }
    }
    
    if (input_filename[0] == 0) {
        printf("Error: you must specify an input filename.\n");
        exit(2);
    }
    if (output_filename[0] == 0) {
        printf("Error: you must specify an output filename.\n");
        exit(2);
    }
    if (database_dir[0] == 0) {
        printf("Error: you must specify a database directory.\n");
        exit(2);
    }
}

/*----------------------------------------------------------------------*
 * Function:   chomp
 * Purpose:    Remove hidden characters from end of line
 * Parameters: str -> String to change
 * Returns:    None
 *----------------------------------------------------------------------*/
void chomp(char* str)
{
    int i = strlen(str) - 1;
    
    while ((i > 0) && (str[i] < ' ')) {
        str[i--] = 0;
    }
}

/*----------------------------------------------------------------------*
 * Function:   allocate_memory
 * Purpose:    Allocate memory to store tabkes
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void allocate_memory(void)
{
    if (is_gi) {
        memory_required+=(max_gi * sizeof(unsigned int));
        printf("Allocating memory for GI list (%d entries)\n", max_gi);
        gi_to_node = calloc(max_gi, sizeof(unsigned int));
        if (!gi_to_node) {
            printf("Error: couldn't allocate memory.\n");
            exit(1);
        }
    }
    
    memory_required+=(MAX_NAMES * sizeof(char*));
    printf("Allocating memory for names list (%d entries)\n", MAX_NAMES);
    names = calloc(MAX_NAMES, sizeof(char*));
    if (!names) {
        printf("Error: couldn't allocate memory.\n");
        exit(1);
    }

    memory_required+=(MAX_NAMES * sizeof(int*));
    printf("Allocating memory for nodes list (%d entries)\n", MAX_NAMES);
    nodes = calloc(MAX_NAMES, sizeof(int*));
    if (!nodes) {
        printf("Error: couldn't allocate memory.\n");
        exit(1);
    }

    printf("Total memory required %ld Mb\n", memory_required / (1024*1025));
}

/*----------------------------------------------------------------------*
 * Function:   load_gi_node_list
 * Purpose:    Load GI to node list translation file
 * Parameters: filename -> filename of file
 * Returns:    None
 *----------------------------------------------------------------------*/
void load_gi_to_node_list()
{
    char filename[MAX_PATH];
    char line[1024];
    FILE* fp;
    
    if (is_nucleotide) {
        sprintf(filename, "%s/gi_taxid_nucl.dmp", database_dir);
    } else {
        sprintf(filename, "%s/gi_taxid_prot.dmp", database_dir);        
    }
    printf("Opening database file %s\n", filename);
    fp = fopen(filename, "r");
    if (!fp) {
        printf("Error: can't open %s\n", filename);
    }    
    
    while (!feof(fp)) {
        if (fgets(line, 1024, fp)) {
            char* gi_str = strtok(line, "\t");
            char* node_id_str = strtok(0, "\t");
            
            if ((!gi_str) || (!node_id_str)) {
                printf("Error: bad line in GI file\n");
            } else {
                unsigned int gi = atoi(gi_str);
                unsigned int node_id = atoi(node_id_str);
                
                if (gi >= max_gi) {
                    printf("Error: GI out of range - %d\n", gi);
                    exit(1);
                } else {
                    gi_to_node[gi] = node_id;
                }
            
            }
        
        }
    }
           
    fclose(fp);
}

/*----------------------------------------------------------------------*
 * Function:   load_node_list
 * Purpose:    Load the list of nodes and parent nodes
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void load_node_list(void)
{
    char filename[MAX_PATH];
    char line[1024];
    FILE* fp;
    
    sprintf(filename, "%s/nodes.dmp", database_dir);
    printf("Opening database file %s\n", filename);
    fp = fopen(filename, "r");
    if (!fp) {
        printf("Error: can't open %s\n", filename);
    }
        
    while (!feof(fp)) {
        if (fgets(line, 1024, fp)) {
            char* child_str = strtok(line, "\t");
            char* ignore = strtok(0, "\t");
            char* parent_str = strtok(0, "\t");
            
            if ((!child_str) || (!parent_str)) {
                printf("Error: bad line in nodes file\n");
            } else {
                unsigned int child = atoi(child_str);
                unsigned int parent = atoi(parent_str);                
                nodes[child] = parent;
            }
        }
    }
    
    fclose(fp);
}

/*----------------------------------------------------------------------*
 * Function:   get_name_fields
 * Purpose:    Get name fields from an entry in the name file. Can't use
 *             simple strtok because sometimes fields are missing.
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void get_name_fields(char* string, char* id, char* name, char* unique_name, char* class)
{
    char fields[16][128];
    int field = 0;
    int p = 0;
    char* s = string;
    
    // Separate fields
    fields[0][0] = 0;    
    while (*s != 0) {
        if (*s == '\t') {
            fields[field][p] = 0;
            field++;
            p = 0;
            fields[field][p] = 0;
        } else {
            fields[field][p++] = *s;
        }
        s++;
    }
    fields[field][p] = 0;

    // Return values
    strcpy(id, fields[0]);
    strcpy(name, fields[2]);
    strcpy(unique_name, fields[4]);
    strcpy(class, fields[6]);
}

/*----------------------------------------------------------------------*
 * Function:   load_name_list
 * Purpose:    Load the name list - converts between node IDs and names
 * Parameters: filename -> filename of name list
 * Returns:    None
 *----------------------------------------------------------------------*/
void load_name_list(void)
{
    char filename[MAX_PATH];
    char line[1024];
    FILE* fp;
    
    sprintf(filename, "%s/names.dmp", database_dir);
    printf("Opening database file %s\n", filename);
    fp = fopen(filename, "r");    
    if (!fp) {
        printf("Error: can't open %s\n", filename);
    }
    
    while (!feof(fp)) {
        if (fgets(line, 1024, fp)) {
            char id_str[128];
            char name[128];
            char unique_name[128];
            char class[128];
            
            get_name_fields(line, id_str, name, unique_name, class);
            
            unsigned int id = atoi(id_str);

            if (strcmp(class, "scientific name") == 0) {
                memory_required += (strlen(name)+1);
                names[id] = malloc(strlen(name)+1);
                if (!names[id]) {
                    printf("Error: can't allocate memory for name\n");
                } else {
                    strcpy(names[id], name);
                }
            }
        }
    }
    
    fclose(fp);
}

/*----------------------------------------------------------------------*
 * Function:   get_taxonomy_by_gi
 * Purpose:    Given a GI, return taxonomy string
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
char* get_taxonomy_from_node(int current_node, char* taxonomy)
{
    int node_list[1024];
    int n = 0;
    int i;

    taxonomy[0] = 0;

    while (current_node > 1) {
        node_list[n++] = current_node;
        current_node = nodes[current_node];
    }
    
    for (i=n-1; i>=0; i--) {
        if (names[node_list[i]] > 0) {
            strcat(taxonomy, names[node_list[i]]);
            if (i > 0) strcat(taxonomy, ",");
        } else {
            strcat(taxonomy, "Unknown");
            if (i > 0) strcat(taxonomy, ",");
            printf("\nError: no name for node %d\n", node_list[i]);
        }
    }
    
    return taxonomy;
}

/*----------------------------------------------------------------------*
 * Function:   get_taxonomy_by_gi
 * Purpose:    Given a GI, return taxonomy string
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
char* get_taxonomy_by_gi(int gi, char* taxonomy)
{
    int e = 0;
    int current_node;
    
    if ((gi >= max_gi) || (gi < 1)) {
        printf("\nError: bad GI (%d)\n", gi);
        e = 1;
    }
    
    if (e == 0) {
        current_node = gi_to_node[gi];
        
        if (current_node == 0) {
            printf("\nError: GI (%d) node (%d) invalid\n", gi, current_node);
            e = 2;
        }
    }
    
    if (e == 0) {
        taxonomy = get_taxonomy_from_node(current_node, taxonomy);
    }
    
    return taxonomy;
}

/*----------------------------------------------------------------------*
 * Function:   process_request_file
 * Purpose:    Find taxonomy for a file of GIs
 * Parameters: filename -> name of file containing one GI per line
 * Returns:    None
 *----------------------------------------------------------------------*/
void process_request_file()
{
    FILE *fp_in;
    FILE *fp_out;
    char buffer[1024];
    char *accession;
    char *version;
    char line[1024];
    char t[1024];
    int count = 0;
    long int taxid;
    long int gi;
    int found;

    fp_in = fopen(input_filename, "r");
    if (!fp_in) {
        printf("Error: can't open %s\n", input_filename);
        exit(1);
    }

    fp_out = fopen(output_filename, "w");
    if (!fp_out) {
        fclose(fp_in);
        printf("Error: can't open %s\n", output_filename);
        exit(1);
    }    
    
    while (!feof(fp_in)) {        
        if (fgets(line, 1024, fp_in)) {
            chomp(line);
            count++;
            if ((count % 100) == 0) {
                printf(".");
                fflush(stdout);
            }

            if (is_gi) {
                int gi = atoi(line);
                
                if (gi < 1) {
                    printf("Error: bad GI (%d) in request file\n", gi);
                } else {
                    get_taxonomy_by_gi(gi, t);
                    fprintf(fp_out, "%i\t%s\n", gi, t);
                }
            } else if (is_accession) {
                char* query = line;
                
                found = find_accession(query, buffer, &accession, &version, &taxid, &gi);
                if (found == 1) {
                    if (taxid == 0) {
                        strcpy(t, "Unknown");
                    } else {
                        get_taxonomy_from_node(taxid, t);
                    }
                    //printf("Found %s: %s, %s, %d, %d\n", query, accession, version, taxid, gi);
                    fprintf(fp_out, "%s\t%s\n", query, t);
                } else {
                    printf("\nCouldn't find: [%s]\n", query);
                }
            }
        }
    }
    
    fclose(fp_out);
    fclose(fp_in);
    
    printf("\n\nDone. Processed %d IDs.\n", count);
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
void open_acc_file(char* filename) {
    acc_fp = fopen(filename, "r");
    if (!acc_fp) {
        printf("Error: can't open %s\n", filename);
        exit(1);
    }
    fseek(acc_fp, 0, SEEK_END);
    acc_file_size = ftell(acc_fp);
    printf("File size: %li\n", acc_file_size);
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
void close_acc_file() {
    if (acc_fp) {
        fclose(acc_fp);
    }
}

char* get_first_token(char* string, char* value, char token) {
    int i;
    
    for (i=0; i<strlen(string); i++) {
        if (string[i] == token) {
            break;
        } else {
            value[i] = string[i];
        }
    }
    
    return value;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
void get_closest_record(long int pos, char* line) {
    char c;
    char temp[1024];
    char token[1024];
    char* this_acc;
    char* prev_acc;
    char* next_acc;
    
    while (pos >= 0) {
        fseek(acc_fp, pos, SEEK_SET);
        pos--;
        
        c = fgetc(acc_fp);
        if (c == '\n') {
            break;
        }
    }
    
    fgets(line, 1024, acc_fp);
    
#ifdef DEBUG
    printf("Got line: %s\n", line);
#endif
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
void split_fields(char* line, char** accession, char** version, long int* taxid, long int* gi) {
    char* taxid_str;
    char* gi_str;
    
    *accession = strtok(line, "\t");
    *version = strtok(NULL, "\t");
    taxid_str = strtok(NULL, "\t");
    gi_str = strtok(NULL, "\t");
    
    if (taxid_str != NULL) {
        *taxid = atoi(taxid_str);
    } else {
        *taxid = 0;
    }
    
    if (gi_str != NULL) {
        *gi = atoi(gi_str);
    } else {
        *gi = 0;
    }
    
#ifdef DEBUG
    printf("Accession: %s\n", *accession);
    printf("Version: %s\n", *version);
    printf("Tax ID: %ld\n", *taxid);
    printf("GI: %ld\n", *gi);
#endif
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
int find_accession(char* search_accession, char* line, char** accession, char** version, long int* taxid, long int* gi) {
    long int min = 0;
    long int max = acc_file_size;
    int similarity;
    int found = 0;
    
#ifdef DEBUG
    printf("Finding %s\n", search_accession);
#endif
    
    while (found == 0) {
        long int current_pos = min + ((max-min) / 2);
#ifdef DEBUG
        printf("Min: %d Max: %d\n", min, max);
#endif
        get_closest_record(current_pos, line);
        split_fields(line, accession, version, taxid, gi);
        
        similarity = strcmp(*accession, search_accession);
        if (similarity == 0) {
            found = 1;
        } else if (similarity > 0) {
            max = current_pos;
        } else if (similarity < 0) {
            min = current_pos;
        }
        
        if (abs(max - min) < 20) {
            break;
        }
    }
    
    return found;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
void load_accession_file(void)
{
    char filename[MAX_PATH];
    
    if (is_nucleotide) {
        sprintf(filename, "%s/acc2tax_nucl_all.txt", database_dir);
    } else {
        sprintf(filename, "%s/acc2tax_prot_all.txt", database_dir);
    }
    
    printf("Opening database file %s\n", filename);
    
    open_acc_file(filename);
}

/*----------------------------------------------------------------------*
 * Function:   main
 *----------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    //setbuf(stdout, NULL);
    
    printf("\nacc2tax %s\n\n", VERSION);

    parse_command_line(argc, argv);
    allocate_memory();
    if (is_gi) {
        load_gi_to_node_list();
    } else if (is_accession) {
        load_accession_file();
    }
    load_node_list();
    load_name_list();

    printf("Memory required: %ld MB\n\n", memory_required / (1024 * 1024));
    
    process_request_file();

    if (is_accession) {
        close_acc_file();
    }
    
    return 0;
}
