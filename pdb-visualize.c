#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <SDL2/SDL.h>
#include <stdbool.h>

typedef struct {
  char name[5];   // Atom name (e.g. "CA")
  char type;      // Atom type (e.g. 'C', 'O', or 'N')
  float radius;
  float x;        // x coordinate
  float y;        // y coordinate
  float z;        // z coordinate
} Atom;


// Function to set the van der Waals radius of an Atom based on its type
void setVDWRadius(Atom* atom) {
  // Set the van der Waals radius based on the atom type
  if (strcmp(atom->name, "H") == 0) {
    atom->radius = 1.2;
  } else if (strcmp(atom->name, "C") == 0) {
    atom->radius = 1.7;
  } else if (strcmp(atom->name, "N") == 0) {
    atom->radius = 1.55;
  } else if (strcmp(atom->name, "O") == 0) {
    atom->radius = 1.52;
  } else if (strcmp(atom->name, "F") == 0) {
    atom->radius = 1.47;
  } else if (strcmp(atom->name, "P") == 0) {
    atom->radius = 1.8;
  } else if (strcmp(atom->name, "S") == 0) {
    atom->radius = 1.8;
  } else if (strcmp(atom->name, "Cl") == 0) {
    atom->radius = 1.75;
  } else if (strcmp(atom->name, "Br") == 0) {
    atom->radius = 1.85;
  } else if (strcmp(atom->name, "I") == 0) {
    atom->radius = 1.98;
  } else {
    // Unknown atom type
    atom->radius = 1.5;
  }
}

// Function to parse a PDB file and store the atom data in a dynamically-allocated array
Atom* parsePDB(const char* filename, int* numAtoms) {
  // Open the PDB file
  FILE* file = fopen(filename, "r");
  if (file == NULL) {
    perror("Error opening file");
    return NULL;
  }

  // Allocate an array to hold the atoms
  Atom* atoms = malloc(sizeof(Atom));
  if (atoms == NULL) {
    perror("Error allocating memory");
    fclose(file);
    return NULL;
  }

  // Read each line of the file
  char line[80];
  while (fgets(line, sizeof(line), file) != NULL) {
    // Check if this is an atom record
    if (strncmp(line, "ATOM", 4) == 0) {
      // Parse the atom data
      Atom atom;
      sscanf(line, "ATOM  %*d %4s %c %*s %*s %f %f %f", atom.name, &atom.type, &atom.x, &atom.y, &atom.z);

      // Add the atom to the array
      atoms = realloc(atoms, sizeof(Atom) * (*numAtoms + 1));
      atoms[*numAtoms] = atom;
      (*numAtoms)++;
    }
  }

  // Close the file and return the atom array
  fclose(file);
  return atoms;
}

void SDL_RenderDrawCircle(SDL_Renderer* renderer, int x, int y, int radius) {
  const int32_t diameter = (radius * 2);

  int32_t x_offset = radius - 1;
  int32_t y_offset = 0;
  int32_t err = 0;

  while (x_offset >= y_offset) {
    SDL_RenderDrawPoint(renderer, x + x_offset, y + y_offset);
    SDL_RenderDrawPoint(renderer, x + y_offset, y + x_offset);
    SDL_RenderDrawPoint(renderer, x - y_offset, y + x_offset);
    SDL_RenderDrawPoint(renderer, x - x_offset, y + y_offset);
    SDL_RenderDrawPoint(renderer, x - x_offset, y - y_offset);
    SDL_RenderDrawPoint(renderer, x - y_offset, y - x_offset);
    SDL_RenderDrawPoint(renderer, x + y_offset, y - x_offset);
    SDL_RenderDrawPoint(renderer, x + x_offset, y - y_offset);

    y_offset++;
    err += 1 + 2 * y_offset;
    if (2 * (err - x_offset) + 1 > 0) {
      x_offset--;
      err += 1 - 2 * x_offset;
    }
  }
}

// Function to draw an atom as a sphere using the SDL library
void drawAtom(SDL_Renderer* renderer, Atom atom, float scale, int offsetX, int offsetY, int offsetZ) {
  // Calculate the screen coordinates of the atom
  int x = (int)(atom.x * scale) + offsetX;
  int y = (int)(atom.y * scale) + offsetY;
  int z = (int)(atom.z * scale) + offsetZ;

  // Set the draw color based on the atom type
  switch (atom.type) {
    case 'H': SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255); break;
    case 'C': SDL_SetRenderDrawColor(renderer, 200, 200, 200, 255); break;
    case 'N': SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255); break;
    case 'O': SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255); break;
    case 'F': SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255); break;
    case 'P': SDL_SetRenderDrawColor(renderer, 255, 165, 0, 255); break;
    case 'S': SDL_SetRenderDrawColor(renderer, 255, 255, 0, 255); break;
    default:
      // Unknown atom type
      SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
      break;
  }

  // Draw the atom as a filled circle
  for (int r = atom.radius * scale; r > 0; r--) {
    SDL_RenderDrawCircle(renderer, x, y, r);
  }
}

void calculateScaleAndOffset(Atom* atoms, int numAtoms, int screenWidth, int screenHeight, int screenDepth, float* scale, int* offsetX, int* offsetY, int* offsetZ) {
  // Find the minimum and maximum coordinates of the atoms
  float minX = atoms[0].x, maxX = atoms[0].x;
  float minY = atoms[0].y, maxY = atoms[0].y;
  float minZ = atoms[0].z, maxZ = atoms[0].z;
  for (int i = 1; i < numAtoms; i++) {
    minX = fmin(minX, atoms[i].x);
    maxX = fmax(maxX, atoms[i].x);
    minY = fmin(minY, atoms[i].y);
    maxY = fmax(maxY, atoms[i].y);
    minZ = fmin(minZ, atoms[i].z);
    maxZ = fmax(maxZ, atoms[i].z);
  }

  // Calculate the size of the molecule in each dimension
  float sizeX = maxX - minX;
  float sizeY = maxY - minY;
  float sizeZ = maxZ - minZ;

  // Calculate the scale factor to fit the molecule on the screen
  float scaleX = screenWidth / sizeX;
  float scaleY = screenHeight / sizeY;
  float scaleZ = screenDepth / sizeZ;
  *scale = fmin(scaleX, fmin(scaleY, scaleZ));

  // Calculate the offset to center the molecule on the screen
  *offsetX = (int)((screenWidth - (sizeX * *scale)) / 2.0f);
  *offsetY = (int)((screenHeight - (sizeY * *scale)) / 2.0f);
  *offsetZ = (int)((screenDepth - (sizeZ * *scale)) / 2.0f);
}


int main(int argc, char* argv[]) {
  // Parse the PDB file and create an array of Atom structures
  Atom* atoms;
  int numAtoms;
  int screenHeight = 1024;
  int screenWidth = 1024;
  int screenDepth = 50;
  atoms = parsePDB(argv[1], &numAtoms);

  // Set the van der Waals radius of each atom based on its type
  for (int i = 0; i < numAtoms; i++) {
    setVDWRadius(&atoms[i]);
  }

  // Initialize SDL
  if (SDL_Init(SDL_INIT_VIDEO) < 0) {
    fprintf(stderr, "SDL initialization failed: %s\n", SDL_GetError());
    return 1;
  }

  // Create the window and renderer
  SDL_Window* window = SDL_CreateWindow("PDB Viewer", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, screenWidth, screenHeight, SDL_WINDOW_SHOWN);
  SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);

  // Set the scale and offsets for the atoms
  float scale = 10.0;
  int offsetX = 500;
  int offsetY = 200;
  int offsetZ = 0;

  // Main loop
  bool running = true;
  while (running) {
    // Handle events
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      switch (event.type) {
        case SDL_QUIT:
          running = false;
          break;
      }
    }

    // Clear the screen
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    // Draw the atoms
    for (int i = 0; i < numAtoms; i++) {
      drawAtom(renderer, atoms[i], scale, offsetX, offsetY, offsetZ);
    }

    // Update the screen
    SDL_RenderPresent(renderer);

    SDL_Delay(10);
  }

  // Clean up
  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  SDL_Quit();

  return 0;
}
