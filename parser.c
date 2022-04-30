#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"runge-kutta-fehlberg.h"
#include"epidemic.h"

void parse_keywords(char* keyword, char* line, model_parameters *mp) {

  size_t offset;
  int adder;
  int i,j;
  
  if (!strncmp(keyword,"runtime",80)) {
    sscanf(line+7,"%lf",&(mp->runtime));
  } else if (!strncmp(keyword,"n_variants",80)) {
    sscanf(line+10,"%li",&(mp->n_variants));
  } else if (!strncmp(keyword,"n_classes",80)) {
    sscanf(line+9,"%li",&(mp->n_classes));
  } else if (!strncmp(keyword,"error_tolerance",80)) {
    sscanf(line+15, "%lf",&(mp->error_tolerance));
  } else if (!strncmp(keyword,"total_population",80)) {
    sscanf(line+16,"%lf",&(mp->total_population));
  } else if (!strncmp(keyword,"population_class_distribution",80)) {
    offset = 29;
    mp->population_class_distribution =
      (double*)malloc(sizeof(double)*mp->n_classes);
    for(i=0;i<mp->n_classes;i++) {
      sscanf(line+offset,"%lf%n",(mp->population_class_distribution)+i,&adder);
      offset += adder;
    }
  } else if (!strncmp(keyword,"contacts_between_classes", 80)) {
    offset = 24;
    mp->contacts_between_classes =
      (double*)malloc(sizeof(double)*mp->n_classes*mp->n_classes);
    for(i=0;i<mp->n_classes;i++) {
      for(j=0;j<mp->n_classes;j++) {
	sscanf(line+offset,"%lf%n",
	       (mp->contacts_between_classes)+(i*mp->n_classes+j),&adder);
	offset+=adder;
      }
    }
  } else if (!strncmp(keyword,"reporting_rate",80)) {
    sscanf(line+14,"%lf",&(mp->reporting_rate));
  } else if (!strncmp(keyword,"epsilon",80)) {
    offset = 7;
    mp->epsilon =
      (double*)malloc(sizeof(double)*mp->n_variants*mp->n_variants);
    for(i=0;i<mp->n_variants;i++) {
      for(j=0;j<mp->n_variants;j++) {
        sscanf(line+offset,"%lf%n",(mp->epsilon)+(i*mp->n_variants+j),&adder);
	offset+=adder;
      }
    }
  } else if (!strncmp(keyword,"recovery_rate",80)) {
    offset = 13;
    mp->recovery_rate =
      (double*)malloc(sizeof(double)*mp->n_variants);
    for(i=0;i<mp->n_variants;i++) {
      sscanf(line+offset,"%lf%n",(mp->recovery_rate)+i,&adder);
      offset+=adder;
    }
  } else if (!strncmp(keyword,"sigma",80)) {
    offset = 5;
    mp->sigma =
      (double*)malloc(sizeof(double)*mp->n_variants*mp->n_classes);
    for(j=0;j<mp->n_variants;j++) {  
      for(i=0;i<mp->n_classes;i++) {
	sscanf(line+offset,"%lf%n",(mp->sigma)+(i*mp->n_variants+j),&adder);
	offset+=adder;
      }
    }
  } else if (!strncmp(keyword,"vaccine_logistic_maximum",80)) {
    offset = 24;
    mp->vaccine_logistic_maximum =
      (double*)malloc(sizeof(double)*mp->n_variants*mp->n_classes);
    for(j=0;j<mp->n_classes;j++) {
      for(i=0;i<mp->n_variants;i++) {
	sscanf(line+offset,"%lf%n",
	       (mp->vaccine_logistic_maximum)+(j*mp->n_variants+i),&adder);
	offset+=adder;
      }
    }
  } else if (!strncmp(keyword,"vaccine_logistic_growth",80)) {
    offset = 23;
    mp->vaccine_logistic_growth =
      (double*)malloc(sizeof(double)*mp->n_variants*mp->n_classes);
    for(j=0;j<mp->n_classes;j++) {
      for(i=0;i<mp->n_variants;i++) {
	sscanf(line+offset,"%lf%n",
	       (mp->vaccine_logistic_growth)+(j*mp->n_variants+i),&adder);
	offset+=adder;
      }
    }
  } else if (!strncmp(keyword,"vaccine_logistic_midpoint",80)) {
    offset = 25;
    mp->vaccine_logistic_midpoint =
      (double*)malloc(sizeof(double)*mp->n_variants*mp->n_classes);
    for(j=0;j<mp->n_classes;j++) {
      for(i=0;i<mp->n_variants;i++) {
	sscanf(line+offset,"%lf%n",
	       (mp->vaccine_logistic_midpoint)+(j*mp->n_variants+i),&adder);
	offset+=adder;
      }
    }
  } else if (!strncmp(keyword,"vaccine_efficiency",80)) {
    offset = 18;
    mp->vaccine_efficiency =
      (double*)malloc(sizeof(double)*mp->n_variants*mp->n_classes);
    for(i=0;i<mp->n_classes;i++) {
      for(j=0;j<mp->n_variants;j++) {  
	sscanf(line+offset,"%lf%n",
	       (mp->vaccine_efficiency)+(i*mp->n_variants+j),&adder);
	offset+=adder;
      }
    }
  } else if(!strncmp(keyword,"conversion_rate",80)) {
    offset = 15;
    mp->conversion_rate =
      (double*)malloc(sizeof(double)*mp->n_variants);
    for(i=0;i<mp->n_variants;i++) {
      sscanf(line+offset,"%lf%n",(mp->conversion_rate+i),&adder);
      offset+=adder;
    }
  } else if(!strncmp(keyword,"asymptomatics_fraction",80)) {
    offset = 22;
    mp->asymptomatics_fraction =
      (double*)malloc(sizeof(double)*mp->n_variants*mp->n_classes);
    for(i=0;i<mp->n_classes;i++) {
      for(j=0;j<mp->n_variants;j++) {
	sscanf(line+offset,"%lf%n",
	       (mp->asymptomatics_fraction)+(i*mp->n_variants+j),&adder);
	offset+=adder;
      }
    }
  } else if(!strncmp(keyword,"initial_intermediate_fraction",80)) {
    offset = 29;
    mp->initial_intermediate_fraction =
      (double*)malloc(sizeof(double)*mp->n_variants*mp->n_classes);
    for(i=0;i<mp->n_classes;i++) {
      for(j=0;j<mp->n_variants;j++) {
	sscanf(line+offset,"%lf%n",
	       (mp->initial_intermediate_fraction)+(i*mp->n_variants+j),&adder);
	offset+=adder;
      }
    }
  } else if(!strncmp(keyword,"initial_infected_fraction",80)) {
    offset = 25;
    mp->initial_infected_fraction =
      (double*)malloc(sizeof(double)*mp->n_variants*mp->n_classes);
    for(i=0;i<mp->n_classes;i++) {
      for(j=0;j<mp->n_variants;j++) {
	sscanf(line+offset,"%lf%n",
	       (mp->initial_infected_fraction)+(i*mp->n_variants+j),&adder);
	offset+=adder;
      }
    }
  } else if(!strncmp(keyword,"pre_immune_fraction",80)) {
    offset = 19;
    mp->pre_immune_fraction =
      (double*)malloc(sizeof(double)*mp->n_variants*mp->n_classes);
    for(i=0;i<mp->n_classes;i++) {
      for(j=0;j<mp->n_variants;j++) {
	sscanf(line+offset,"%lf%n",
	       (mp->pre_immune_fraction)+(i*mp->n_variants+j),&adder);
	offset+=adder;
      }
    }

  } else if (!strncmp(keyword,"Rzero",80)) {
    offset = 5;
    mp->Rzero =
      (double*)malloc(sizeof(double)*mp->n_variants);
    for(i=0;i<mp->n_variants;i++) {
      sscanf(line+offset,"%lf%n",(mp->Rzero)+i,&adder);
      offset+=adder;
    }
  } else if (!strncmp(keyword,"variant_introduction_date",80)) {
    offset = 25;
    mp->variant_introduction_date =
      (double*)malloc(sizeof(double)*mp->n_variants);
    for(i=0;i<mp->n_variants;i++) {
      sscanf(line+offset,"%lf%n",(mp->variant_introduction_date)+i,&adder);
      offset+=adder;
    }
  } else if (!strncmp(keyword,"immune_drain",80)) {
    offset = 12;
    mp->immune_drain =
      (double*)malloc(sizeof(double)*mp->n_variants*mp->n_classes);
    for(i=0;i<mp->n_classes;i++) {
      for(j=0;j<mp->n_variants;j++) {  
	sscanf(line+offset,"%lf%n",
	       (mp->immune_drain)+(i*mp->n_variants+j),
	       &adder);
	offset+=adder;
      }
    }
  } else if (!strncmp(keyword,"vaccine_drain",80)) {
    offset = 13;
    mp->vaccine_drain =
      (double*)malloc(sizeof(double)*mp->n_variants*mp->n_classes);
    for(i=0;i<mp->n_classes;i++) {
	for(j=0;j<mp->n_variants;j++) {  

	sscanf(line+offset,"%lf%n",
	       (mp->vaccine_drain)+(i*mp->n_variants+j),
	       &adder);
	offset+=adder;
      }
    }
  }

} 

void parse_model_properties(FILE* in, model_parameters* mp) {

  char* line = NULL;
  size_t n = 0;

  char* keyword = (char*)malloc(sizeof(char)*80);
  
  while ( -1 != getline(&line,&n,in) ) {
    if(1 == sscanf(line,"%s",keyword)) {
      parse_keywords(keyword, line, mp);
    }
  }
  free(keyword);
}
  
    
