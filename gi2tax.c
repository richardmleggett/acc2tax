/*----------------------------------------------------------------------*
 * File:    gi2tax.c
 * Author:  Richard Leggett (richard.leggett@tgac.ac.uk)
 * Purpose: Convert GI number to taxonomy
 * Created: 10 Jan 2013
 *----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h> 

/*----------------------------------------------------------------------*
 * Constants
 *----------------------------------------------------------------------*/
#define MAX_GI 500000000
#define MAX_NAMES 2000000
#define MAX_PATH 10000

/*----------------------------------------------------------------------*
 * Globals
 *----------------------------------------------------------------------*/
char database_dir[MAX_PATH];
char input_filename[MAX_PATH];
char output_filename[MAX_PATH];
int is_nucleotide = 1;
int is_protein = 0;
long int memory_required = 0;
unsigned int* gi_to_node;
char** names;
unsigned int* nodes;

/*----------------------------------------------------------------------*
 * Function:   usage
 * Purpose:    Report program usage.
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void usage(void)
{
    printf("\ngi2tax v0.1\n" \
           "richard.leggett@tgac.ac.uk\n" \
           "\nProvide batch taxonomy information for Genbank iDs.\n" \
           "\nOptions:\n" \
           "    [-h | --help] This help screen.\n" \
           "    [-i | --input] File of Genbank IDs, one per line." \
           "    [-o | --output] Filename of output file." \
           "    [-d | --database] Directory containing NCBI taxonomy files." \
           "    [-n | --nucleotide] Query GIs are nucleotide [default]." \
           "    [-p | --protein] Query GIs are protein." \
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
        {"database", required_argument, NULL, 'd'},
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
    
    while ((opt = getopt_long(argc, argv, "d:hi:no:p", long_options, &longopt_index)) > 0)
    {
        switch(opt) {
            case 'h':
                usage();
                exit(0);
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
 * Function:   allocate_memory
 * Purpose:    Allocate memory to store tabkes
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void allocate_memory(void)
{
    printf("Allocating memory for GI list\n");
    memory_required+=(MAX_GI * sizeof(unsigned int));
    gi_to_node = calloc(MAX_GI, sizeof(int));
    if (!gi_to_node) {
        printf("Error: couldn't allocate memory.\n");
        exit(1);
    }
    
    printf("Allocating memory for names list\n");
    memory_required+=(MAX_NAMES * sizeof(char*));
    names = calloc(MAX_NAMES, sizeof(char*));
    if (!names) {
        printf("Error: couldn't allocate memory.\n");
        exit(1);
    }

    printf("Allocating memory for nodes list\n");
    memory_required+=(MAX_NAMES * sizeof(int*));
    nodes = calloc(MAX_NAMES, sizeof(int*));
    if (!nodes) {
        printf("Error: couldn't allocate memory.\n");
        exit(1);
    }
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
                
                if ((gi < 0) || (gi >= MAX_GI)) {
                    printf("Error: GI out of range - %d\n", gi);
                } else {
                    if (node_id < 0) {
                        printf("Error: Node ID out of range - %d\n", node_id);
                    } else {
                        gi_to_node[gi] = node_id;
                    }
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
char* get_taxonomy_by_gi(int gi, char* taxonomy)
{
    int e = 0;
    int node_list[1024];
    int current_node;
    int n = 0;
    int i;
    
    taxonomy[0] = 0;
    
    if ((gi >= MAX_GI) || (gi < 1)) {
        printf("Error: bad GI (%d)\n", gi);
        e = 1;
    }
    
    if (e == 0) {
        current_node = gi_to_node[gi];
        
        if (current_node == 0) {
            printf("Error: GI (%d) node (%d) invalid\n", gi, current_node);
            e = 2;
        }
    }
    
    while ((current_node > 1) && (e == 0)) {
        node_list[n++] = current_node;
        current_node = nodes[current_node];
    }
    
    for (i=n-1; i>=0; i--) {
        if (names[node_list[i]] > 0) {
            strcat(taxonomy, names[node_list[i]]);
            if (i > 0) strcat(taxonomy, ",");
        } else {
            printf("Error: no name for node %d\n", node_list[i]);
        }
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
    char line[1024];
    char t[1024];
    int count = 0;

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
            int gi = atoi(line);
            
            count++;
            
            if (gi < 1) {
                printf("Error: bad GI (%s) in request file\n", gi);
            } else {
                get_taxonomy_by_gi(gi, t);
                fprintf(fp_out, "%i\t%s\n", gi, t);
            }
        }
    }
    
    fclose(fp_out);
    fclose(fp_in);
    
    printf("Processed %d IDs.\n", count);
}

/*----------------------------------------------------------------------*
 * Function:   main
 *----------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    parse_command_line(argc, argv);
    allocate_memory();
    load_gi_to_node_list();
    load_node_list();
    load_name_list();

    printf("Memory required: %d MB\n\n", memory_required / (1024 * 1024));
    
    process_request_file();
}
