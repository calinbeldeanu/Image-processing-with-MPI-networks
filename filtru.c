#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void topologyFromFile(FILE *topology ,int rank,int *parent, int *NrNeigh, int *neigh) {
	
	*NrNeigh = 0;
	*parent = -1;

	int node = -1;
	char buffer[512];
	char delim[3] = ": ";
	char *token;

	while (fgets (buffer, 512, topology) != NULL) {

    	token = strtok(buffer, delim);
    	sscanf(token, "%d", &node);
    	if (rank == node) {

    		token = strtok(NULL, delim);
			while (token != NULL) {
				sscanf(token , "%d" , &node);
				neigh[*NrNeigh] = node;
				*NrNeigh = *NrNeigh + 1;

				token = strtok(NULL, delim);
			}

			break;
    	}
 	}

 	
 	if(rank != 0){
 		if(*NrNeigh == 1){
 			*parent = neigh[0];
 			*NrNeigh = 0;
 		}else{
 			*parent = neigh[0];
 			for(int i = 0 ; i < *NrNeigh -1 ; i++){
 				neigh[i] = neigh[i+1];
 			}
 			*NrNeigh -=1;
 		}	
 	}
	
 	fclose(topology);
	MPI_Barrier (MPI_COMM_WORLD);
}	

void applyFilter(char **filterType, char **picName, char **newPicName, int rank,int size, int *neigh,int *statistics, int NrNeigh, int parent) {

	FILE *pic, *newPic;
	MPI_Status stat;
	int one = 1;
	int equalpart, reminder, lastpart ,val; 
	int **picMatrix, **newpicMatrix;
	int i , j , l ,k ,m ,start= 0;
	char *token = malloc(128);
	int w, h ,isSobel;
	char buffer[512];
	char buffer1[512];

	for (l = 0; l < size; ++l) {

		if (rank == 0) {
						
			pic = fopen(picName[l], "r");

				fgets (buffer, 512, pic);

				fgets (buffer, 512, pic);
				fgets (buffer, 512, pic);

				token = strtok(buffer, " \n"); 
				sscanf(token , "%d" , &w);
				token = strtok(NULL, " \n"); 
				sscanf(token , "%d" , &h); 

				int *v = (int *) calloc((h+2)*(w+2), sizeof(int)); ;
				int **picMatrix = (int **) calloc((h+2), sizeof(int *));
				for (i = 0; i < h+2; i++) {
						picMatrix[i] = v + i * (w +2);
				}
			
				fgets (buffer, sizeof(buffer), pic);

				for (i = 1; i < h + 1; ++i) {
					for (j = 1; j < w + 1; ++j) {
						fgets (buffer1, 512, pic);

						sscanf(buffer1 , "%d" , &val);
						picMatrix[i][j] = val;	
					}
				}
				fclose(pic);
			
			if (strncmp(filterType[l], "sobel", 4) == 0) {
				isSobel = 1;
			}else{
				isSobel = 0;
			}
			
			equalpart,lastpart = 0;
			equalpart = h / NrNeigh;
			reminder = h % NrNeigh;
			lastpart = equalpart + reminder;

			int *v1 = (int *) calloc((h)*(w), sizeof(int)); ;
			int **newpicMatrix = (int **) calloc((h)*1, sizeof(int *));
			for (i = 0; i < h; i++) {
				newpicMatrix[i] = v1 + (i * w);
			}

			for(i = 0; i < NrNeigh; i++)
			{	
				if(i != NrNeigh -1){
					MPI_Send(&equalpart, one, MPI_INT, neigh[i], rank, MPI_COMM_WORLD); 
					MPI_Send(&w, one, MPI_INT, neigh[i], rank, MPI_COMM_WORLD); 
					MPI_Send(&isSobel, one, MPI_INT, neigh[i], rank, MPI_COMM_WORLD); 
					MPI_Send(&picMatrix[i*equalpart][start], (equalpart + 2) * (w + 2), MPI_INT, neigh[i], rank, MPI_COMM_WORLD); 
				}else {
					MPI_Send(&lastpart, one, MPI_INT, neigh[i], rank, MPI_COMM_WORLD); 
					MPI_Send(&w, one, MPI_INT, neigh[i], rank, MPI_COMM_WORLD); 
					MPI_Send(&isSobel, one, MPI_INT, neigh[i], rank, MPI_COMM_WORLD); 
					MPI_Send(&picMatrix[i*equalpart][start], (lastpart + 2) * (w + 2), MPI_INT, neigh[i], rank, MPI_COMM_WORLD);
				}
			}
 
		

			for(i = 0; i < NrNeigh ; i++) {
				if(i != NrNeigh -1){
					MPI_Recv(&newpicMatrix[i*equalpart][0], equalpart * w, MPI_INT, neigh[i], MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
				}else{
					MPI_Recv(&newpicMatrix[(i)*equalpart][0], lastpart * w, MPI_INT, neigh[i], MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

				}
			}

			for (i = 0; i < h; i++){
				if(i == 0){
					newPic = fopen(newPicName[l], "w");
					fprintf(newPic, "P2\n");
					fprintf(newPic, "# CREATOR: GIMP PNM Filter Version 1.1\n");
					fprintf(newPic, "%d %d\n", w, h);
					fprintf(newPic, "255\n");
				}
				for (j = 0; j < w; j++){
					fprintf(newPic, "%d\n", newpicMatrix[i][j]);
				}
			}	

			fclose(newPic);

		}else{

			MPI_Recv(&h, one, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD , &stat); 
			MPI_Recv(&w, one, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD , &stat);

			int *v = (int *) calloc((h+2)*(w+2), sizeof(int)); ;
			int **picMatrix = (int **) calloc((h+2), sizeof(int *));
			for (i = 0; i < h+2; i++) {
				picMatrix[i] = v + i * (w +2);
			}

			MPI_Recv(&isSobel, one, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat); 
			MPI_Recv(&picMatrix[start][start], (h + 2) * (w + 2), MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat); 

			equalpart ,lastpart = 0;
			
			if (NrNeigh != 0) {

				int *v1 = (int *) calloc((h*w), sizeof(int)); ;
				int **newpicMatrix = (int **) calloc((h*1), sizeof(int *));
				for (i = 0; i < h; i++) {
					newpicMatrix[i] = v1 + (i * w);
				}

				
				equalpart = h / NrNeigh;
				reminder = h % NrNeigh;
				lastpart = equalpart + reminder;



				for(i = 0; i < NrNeigh; i++)
				{	
					if(i != NrNeigh -1){
						MPI_Send(&equalpart, one, MPI_INT, neigh[i], rank, MPI_COMM_WORLD);
						MPI_Send(&w, one, MPI_INT, neigh[i], rank, MPI_COMM_WORLD);
						MPI_Send(&isSobel, one, MPI_INT, neigh[i], rank, MPI_COMM_WORLD);
						MPI_Send(&picMatrix[i*equalpart][0], (equalpart + 2) * (w + 2), MPI_INT, neigh[i], rank, MPI_COMM_WORLD);
					}else{
						MPI_Send(&lastpart, one, MPI_INT, neigh[i], rank, MPI_COMM_WORLD);
						MPI_Send(&w, one, MPI_INT, neigh[i], rank, MPI_COMM_WORLD);
						MPI_Send(&isSobel, one, MPI_INT, neigh[i], rank, MPI_COMM_WORLD);
						MPI_Send(&picMatrix[i*equalpart][start], (lastpart + 2) * (w + 2), MPI_INT, neigh[i], rank, MPI_COMM_WORLD);

					}
				}
	
				for(i = 0; i < NrNeigh ; i++) {
					if(i != NrNeigh){
						MPI_Recv(&newpicMatrix[i*equalpart][start], equalpart * w, MPI_INT, neigh[i], MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
					}else{
						MPI_Recv(&newpicMatrix[(i)*equalpart][start], lastpart * w, MPI_INT, neigh[i], MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
						equalpart, lastpart = 0;
					}
				}

				MPI_Send(&newpicMatrix[start][start], h*w, MPI_INT, parent, rank, MPI_COMM_WORLD);//atentie

			} else {
				int add = 127;
				int help[9] = {1 , 0,-1,2,0,-2,1,0,-1};
				int help1[9] = {-1,-1,-1,-1,9,-1,-1,-1,-1};
				*statistics = *statistics + h;
			

				int *v1 = (int *) calloc((h*w), sizeof(int)); ;
				int **newpicMatrix = (int **) calloc(h*1, sizeof(int *));
				for (i = 0; i < h; i++) {
					newpicMatrix[i] = v1 + (i * w);
				}

				if(isSobel == 1) {
					i= 1;
					while ( i < h+1){
						for (j = 1; j < w+1; j++){


							newpicMatrix[i-1][j-1] =  
													  help[0] * picMatrix[i-1][j-1] +
													  help[1] * picMatrix[i-1][j] +
													  help[2] * picMatrix[i-1][j+1] + help[3] * picMatrix[i][j-1] +
													  help[4] * picMatrix[i][j] +     help[5] * picMatrix[i][j+1] +
													  help[6] * picMatrix[i+1][j-1] +  help[7] * picMatrix[i+1][j] +
													  help[8] * picMatrix[i+1][j+1];
							
							newpicMatrix[i-1][j-1] += add;


							if (newpicMatrix[i-1][j-1] < 0) {
								newpicMatrix[i-1][j-1] = 0;
							} 

							if (newpicMatrix[i-1][j-1] > 255) {
								newpicMatrix[i-1][j-1] = 255;
							}
						}
						i++;	
					}
				}
				else {
					i = 1;
					while(i < h+1){
						for (j = 1; j < w+1; j++){

							
							newpicMatrix[i-1][j-1] = 
													 help1[0] * picMatrix[i-1][j-1] + 
													 help1[1] * picMatrix[i-1][j] +
													 help1[2] * picMatrix[i-1][j+1] + help1[3] * picMatrix[i][j-1] +
													 help1[4] * picMatrix[i][j] +     help1[5] * picMatrix[i][j+1] +
													 help1[6] * picMatrix[i+1][j-1] + help1[7] * picMatrix[i+1][j] +
													 help1[8] * picMatrix[i+1][j+1];

							
							if (newpicMatrix[i-1][j-1] < 0) {
								newpicMatrix[i-1][j-1] = 0;
							}

							if (newpicMatrix[i-1][j-1] > 255) {
								newpicMatrix[i-1][j-1] = 255;
							}
						}
						i++;	
					}
				}
				
				MPI_Send(&newpicMatrix[0][0], h*w, MPI_INT, parent, rank, MPI_COMM_WORLD); //atentie
			}
		}
		MPI_Barrier (MPI_COMM_WORLD);
	}
	MPI_Barrier (MPI_COMM_WORLD);
}


int main(int argc, char * argv[]) {

	MPI_Init(&argc, &argv);
	MPI_Request request;
	MPI_Status stat;

	int nProcesses,rank , size , parent ,NrNeigh;
	char **filterType ,**picName , **newPicName;
	int *neigh = (int *) calloc(nProcesses, sizeof(int));
	int statistics = 0;
	int *finalvec = (int *) calloc(nProcesses, sizeof(int));
	int *tokens = (int *) calloc(nProcesses, sizeof(int));
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

	FILE *topology = fopen(argv[1], "r");
	topologyFromFile(topology ,rank ,&parent,&NrNeigh,neigh);

	

	if (rank == 0) {
		FILE *images = fopen(argv[2], "r"); 
		char buffer[512];
		char *token = malloc(100);

		fgets (buffer, 512, images);
		sscanf(buffer ,"%d" , &size);

		filterType = (char **) calloc(size, sizeof(char *));
		picName = (char **) calloc(size, sizeof(char *));
		newPicName = (char **) calloc(size, sizeof(char *));
		for (int i = 0; i < size; i++) {
			filterType[i] = (char *) calloc(128, sizeof(char));
			picName[i] = (char *) calloc(128, sizeof(char));
			newPicName[i] = (char *) calloc(128, sizeof(char));
		}

		for (int i = 0; i < size; ++i) {
			fgets (buffer, 512, images);

			token = strtok(buffer, " \n");
			sscanf(token, "%s" , filterType[i]);

			token = strtok(NULL, " \n");
			sscanf(token, "%s" , picName[i]);

			token = strtok(NULL, " \n");	
			sscanf(token, "%s" , newPicName[i]);
		}
		fclose(images);

	} else {
		FILE *images = fopen(argv[2], "r"); 
		char buffer[512];
		char *token = malloc(100);

		fgets (buffer, 512, images);
		sscanf(buffer ,"%d" , &size);

		fclose(images);
	}

	MPI_Barrier (MPI_COMM_WORLD);
	applyFilter(filterType, picName,newPicName,rank,size, neigh, &statistics ,NrNeigh, parent);

	


	if (rank == 0) {

		int j = 0;
		while( j < NrNeigh) {

			MPI_Recv(&tokens[0], nProcesses, MPI_INT, neigh[j], MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

			for(int m = 0; m < nProcesses; m++) {
				finalvec[m] += tokens[m];
			}
			j++;
		}


		FILE *stati = fopen(argv[3], "w"); 
		int p = 0 ; 
		while(p < nProcesses){
			fprintf(stati, "%d: %d\n", p, finalvec[p]);
			p++;
		}

	fclose(stati);

	}else{
 
		if (NrNeigh != 0) {
			
			int i = 0;
			while(i < NrNeigh) {
				MPI_Recv(&tokens[0], nProcesses, MPI_INT, neigh[i], MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

				for(int k = 0; k < nProcesses; k++) {
					finalvec[k] += tokens[k];
				}
				i++;
			}
		
			MPI_Send(&finalvec[0], nProcesses, MPI_INT, parent, rank, MPI_COMM_WORLD); 

		} else {

			finalvec[rank] = statistics ;

			MPI_Send(&finalvec[0], nProcesses, MPI_INT, parent, rank, MPI_COMM_WORLD); 
		}
	}

	MPI_Barrier (MPI_COMM_WORLD);



	MPI_Finalize();
	return 0;
}
