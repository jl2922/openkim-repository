   switch(GetComputeIndex(isComputeProcess_dEdr,
                          isComputeProcess_d2Edr2,
                          isComputeEnergy,
                          isComputeForces,
                          isComputeParticleEnergy,
                          isShift))
   {
      case 0:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 2:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 3:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 4:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 5:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 6:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 7:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 8:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 9:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 10:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 11:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 12:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 13:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 14:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 15:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 16:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 17:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 18:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 19:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 20:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 21:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 22:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 23:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 24:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 25:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 26:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 27:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 28:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 29:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 30:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 31:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 32:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 33:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 34:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 35:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 36:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 37:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 38:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 39:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 40:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 41:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 42:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 43:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 44:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 45:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 46:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 47:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 48:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 49:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 50:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 51:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 52:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 53:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 54:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 55:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 56:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 57:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 58:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 59:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 60:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 61:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 62:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 63:
         ier = Compute< ClusterIterator, false, Coordinates,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 64:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 65:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 66:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 67:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 68:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 69:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 70:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 71:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 72:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 73:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 74:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 75:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 76:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 77:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 78:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 79:
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 80:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 81:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 82:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 83:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 84:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 85:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 86:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 87:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 88:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 89:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 90:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 91:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 92:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 93:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 94:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 95:
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 96:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 97:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 98:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 99:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 100:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 101:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 102:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 103:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 104:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 105:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 106:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 107:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 108:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 109:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 110:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 111:
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 112:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 113:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 114:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 115:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 116:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 117:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 118:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 119:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 120:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 121:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 122:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 123:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 124:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 125:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 126:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 127:
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 128:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 129:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 130:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 131:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 132:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 133:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 134:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 135:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 136:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 137:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 138:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 139:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 140:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 141:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 142:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 143:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 144:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 145:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 146:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 147:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 148:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 149:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 150:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 151:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 152:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 153:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 154:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 155:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 156:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 157:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 158:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 159:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 160:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 161:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 162:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 163:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 164:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 165:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 166:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 167:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 168:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 169:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 170:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 171:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 172:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 173:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 174:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 175:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 176:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 177:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 178:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 179:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 180:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 181:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 182:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 183:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 184:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 185:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 186:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 187:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 188:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 189:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 190:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 191:
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 192:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 193:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 194:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 195:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 196:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 197:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 198:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 199:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 200:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 201:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 202:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 203:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 204:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 205:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 206:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 207:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 208:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 209:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 210:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 211:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 212:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 213:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 214:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 215:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 216:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 217:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 218:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 219:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 220:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 221:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 222:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 223:
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 224:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 225:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 226:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 227:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 228:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 229:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 230:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 231:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 232:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 233:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 234:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 235:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 236:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 237:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 238:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 239:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 240:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 241:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 242:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 243:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 244:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 245:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 246:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 247:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 248:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 249:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 250:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 251:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 252:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 253:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 254:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 255:
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 256:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 257:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 258:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 259:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 260:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 261:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 262:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 263:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 264:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 265:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 266:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 267:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 268:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 269:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 270:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 271:
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 272:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 273:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 274:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 275:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 276:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 277:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 278:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 279:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 280:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 281:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 282:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 283:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 284:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 285:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 286:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 287:
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 288:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 289:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 290:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 291:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 292:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 293:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 294:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 295:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 296:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 297:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 298:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 299:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 300:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 301:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 302:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 303:
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 304:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 305:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 306:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 307:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 308:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 309:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 310:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 311:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 312:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 313:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 314:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 315:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 316:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 317:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 318:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 319:
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 320:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 321:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 322:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 323:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 324:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 325:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 326:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 327:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 328:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 329:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 330:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 331:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 332:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 333:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 334:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 335:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 336:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 337:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 338:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 339:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 340:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 341:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 342:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 343:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 344:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 345:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 346:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 347:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 348:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 349:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 350:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 351:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 352:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 353:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 354:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 355:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 356:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 357:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 358:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 359:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 360:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 361:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 362:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 363:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 364:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 365:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 366:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 367:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 368:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 369:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 370:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 371:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 372:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 373:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 374:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 375:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 376:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 377:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 378:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 379:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 380:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 381:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 382:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 383:
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 384:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 385:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 386:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 387:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 388:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 389:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 390:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 391:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 392:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 393:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 394:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 395:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 396:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 397:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 398:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 399:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 400:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 401:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 402:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 403:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 404:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 405:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 406:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 407:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 408:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 409:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 410:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 411:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 412:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 413:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 414:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 415:
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 416:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 417:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 418:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 419:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 420:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 421:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 422:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 423:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 424:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 425:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 426:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 427:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 428:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 429:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 430:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 431:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 432:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 433:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 434:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 435:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 436:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 437:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 438:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 439:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 440:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 441:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 442:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 443:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 444:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 445:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 446:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 447:
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 448:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 449:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 450:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 451:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 452:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 453:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 454:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 455:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 456:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 457:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 458:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 459:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 460:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 461:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 462:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 463:
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 464:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 465:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 466:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 467:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 468:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 469:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 470:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 471:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 472:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 473:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 474:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 475:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 476:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 477:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 478:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 479:
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 480:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 481:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 482:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 483:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 484:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 485:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 486:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 487:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 488:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 489:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 490:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 491:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 492:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 493:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 494:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 495:
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 496:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 497:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 498:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 499:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 500:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 501:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 502:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 503:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 504:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 505:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 506:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 507:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 508:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 509:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 510:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 511:
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 512:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 513:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 514:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 515:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 516:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 517:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 518:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 519:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 520:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 521:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 522:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 523:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 524:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 525:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 526:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 527:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 528:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 529:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 530:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 531:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 532:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 533:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 534:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 535:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 536:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 537:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 538:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 539:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 540:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 541:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 542:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 543:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 544:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 545:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 546:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 547:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 548:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 549:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 550:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 551:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 552:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 553:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 554:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 555:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 556:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 557:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 558:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 559:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 560:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 561:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 562:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 563:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 564:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 565:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 566:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 567:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 568:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 569:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 570:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 571:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 572:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 573:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 574:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 575:
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 576:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 577:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 578:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 579:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 580:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 581:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 582:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 583:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 584:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 585:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 586:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 587:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 588:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 589:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 590:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 591:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 592:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 593:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 594:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 595:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 596:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 597:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 598:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 599:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 600:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 601:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 602:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 603:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 604:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 605:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 606:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 607:
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 608:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 609:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 610:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 611:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 612:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 613:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 614:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 615:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 616:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 617:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 618:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 619:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 620:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 621:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 622:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 623:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 624:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 625:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 626:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 627:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 628:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 629:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 630:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 631:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 632:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 633:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 634:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 635:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 636:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 637:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 638:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 639:
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 640:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 641:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 642:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 643:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 644:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 645:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 646:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 647:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 648:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 649:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 650:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 651:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 652:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 653:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 654:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 655:
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 656:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 657:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 658:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 659:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 660:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 661:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 662:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 663:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 664:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 665:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 666:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 667:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 668:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 669:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 670:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 671:
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 672:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 673:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 674:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 675:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 676:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 677:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 678:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 679:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 680:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 681:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 682:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 683:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 684:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 685:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 686:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 687:
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 688:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 689:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 690:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 691:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 692:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 693:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 694:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 695:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 696:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 697:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 698:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 699:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 700:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 701:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 702:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 703:
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 704:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 705:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 706:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 707:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 708:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 709:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 710:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 711:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 712:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 713:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 714:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 715:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 716:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 717:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 718:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 719:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 720:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 721:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 722:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 723:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 724:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 725:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 726:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 727:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 728:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 729:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 730:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 731:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 732:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 733:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 734:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 735:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 736:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 737:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 738:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 739:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 740:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 741:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 742:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 743:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 744:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 745:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 746:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 747:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 748:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 749:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 750:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 751:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 752:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 753:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 754:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 755:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 756:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 757:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 758:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 759:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 760:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 761:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 762:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 763:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 764:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 765:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 766:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 767:
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 768:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 769:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 770:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 771:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 772:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 773:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 774:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 775:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 776:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 777:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 778:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 779:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 780:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 781:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 782:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 783:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 784:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 785:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 786:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 787:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 788:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 789:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 790:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 791:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 792:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 793:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 794:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 795:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 796:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 797:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 798:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 799:
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 800:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 801:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 802:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 803:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 804:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 805:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 806:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 807:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 808:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 809:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 810:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 811:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 812:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 813:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 814:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 815:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 816:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 817:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 818:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 819:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 820:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 821:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 822:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 823:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 824:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 825:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 826:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 827:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 828:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 829:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 830:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 831:
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 832:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 833:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 834:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 835:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 836:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 837:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 838:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 839:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 840:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 841:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 842:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 843:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 844:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 845:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 846:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 847:
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 848:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 849:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 850:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 851:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 852:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 853:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 854:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 855:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 856:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 857:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 858:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 859:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 860:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 861:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 862:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 863:
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 864:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 865:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 866:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 867:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 868:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 869:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 870:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 871:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 872:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 873:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 874:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 875:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 876:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 877:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 878:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 879:
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 880:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 881:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 882:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 883:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 884:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 885:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 886:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 887:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 888:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 889:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 890:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 891:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 892:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 893:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 894:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 895:
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 896:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 897:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 898:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 899:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 900:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 901:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 902:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 903:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 904:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 905:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 906:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 907:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 908:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 909:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 910:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 911:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 912:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 913:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 914:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 915:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 916:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 917:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 918:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 919:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 920:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 921:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 922:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 923:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 924:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 925:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 926:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 927:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 928:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 929:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 930:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 931:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 932:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 933:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 934:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 935:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 936:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 937:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 938:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 939:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 940:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 941:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 942:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 943:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 944:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 945:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 946:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 947:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 948:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 949:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 950:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 951:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 952:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 953:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 954:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 955:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 956:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 957:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 958:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 959:
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 960:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 961:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 962:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 963:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 964:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 965:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 966:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 967:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 968:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 969:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 970:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 971:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 972:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 973:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 974:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 975:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 976:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 977:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 978:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 979:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 980:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 981:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 982:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 983:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 984:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 985:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 986:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 987:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 988:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 989:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 990:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 991:
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 992:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 993:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 994:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 995:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 996:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 997:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 998:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 999:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1000:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1001:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1002:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1003:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1004:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1005:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1006:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1007:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1008:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1009:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1010:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1011:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1012:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1013:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1014:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1015:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1016:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1017:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1018:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1019:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1020:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1021:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1022:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1023:
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1024:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1025:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1026:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1027:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1028:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1029:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1030:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1031:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1032:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1033:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1034:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1035:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1036:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1037:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1038:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1039:
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1040:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1041:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1042:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1043:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1044:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1045:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1046:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1047:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1048:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1049:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1050:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1051:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1052:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1053:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1054:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1055:
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1056:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1057:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1058:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1059:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1060:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1061:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1062:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1063:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1064:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1065:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1066:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1067:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1068:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1069:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1070:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1071:
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1072:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1073:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1074:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1075:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1076:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1077:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1078:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1079:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1080:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1081:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1082:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1083:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1084:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1085:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1086:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1087:
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1088:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1089:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1090:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1091:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1092:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1093:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1094:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1095:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1096:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1097:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1098:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1099:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1100:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1101:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1102:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1103:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1104:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1105:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1106:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1107:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1108:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1109:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1110:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1111:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1112:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1113:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1114:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1115:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1116:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1117:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1118:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1119:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1120:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1121:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1122:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1123:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1124:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1125:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1126:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1127:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1128:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1129:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1130:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1131:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1132:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1133:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1134:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1135:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1136:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        false, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1137:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        false, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1138:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        false, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1139:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        false, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1140:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        false, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1141:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        false, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1142:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        false, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1143:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        false, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1144:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        true, false,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1145:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        true, false,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1146:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        true, false,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1147:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        true, false,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1148:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        true, true,
                        false, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1149:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        true, true,
                        false, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1150:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        true, true,
                        true, false >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      case 1151:
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        true, true,
                        true, true >(
                  pkim,
                  particleSpecies,
                  get_neigh,
                  boxSideLengths,
                  coordinates,
                  energy,
                  forces,
                  particleEnergy);
         break;
      default:
         std::cout << "Unknown compute function index" << std::endl;
         ier = KIM_STATUS_FAIL;
         break;
   }
