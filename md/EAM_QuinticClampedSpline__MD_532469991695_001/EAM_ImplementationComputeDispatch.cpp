   switch(GetComputeIndex(isComputeProcess_dEdr,
                          isComputeProcess_d2Edr2,
                          isComputeEnergy,
                          isComputeForces,
                          isComputeParticleEnergy))
   {
      case 0:
         ier = Compute< ClusterIterator, false, Coordinates,
                        false, false,
                        false, false, false >(
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
                        false, false, true >(
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
                        false, true, false >(
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
                        false, true, true >(
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
                        true, false, false >(
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
                        true, false, true >(
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
                        true, true, false >(
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
                        true, true, true >(
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
                        false, true,
                        false, false, false >(
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
                        false, true,
                        false, false, true >(
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
                        false, true,
                        false, true, false >(
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
                        false, true,
                        false, true, true >(
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
                        false, true,
                        true, false, false >(
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
                        false, true,
                        true, false, true >(
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
                        false, true,
                        true, true, false >(
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
                        false, true,
                        true, true, true >(
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
                        true, false,
                        false, false, false >(
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
                        true, false,
                        false, false, true >(
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
                        true, false,
                        false, true, false >(
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
                        true, false,
                        false, true, true >(
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
                        true, false,
                        true, false, false >(
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
                        true, false,
                        true, false, true >(
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
                        true, false,
                        true, true, false >(
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
                        true, false,
                        true, true, true >(
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
                        true, true,
                        false, false, false >(
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
                        true, true,
                        false, false, true >(
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
                        true, true,
                        false, true, false >(
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
                        true, true,
                        false, true, true >(
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
                        true, true,
                        true, false, false >(
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
                        true, true,
                        true, false, true >(
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
                        true, true,
                        true, true, false >(
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
                        true, true,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, false,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        false, true,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, false,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, false, RVec,
                        true, true,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, false,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        false, true,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, false,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, false, MI_OPBC,
                        true, true,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, false,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        false, true,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, false,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, true, Coordinates,
                        true, true,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, false,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        false, true,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, false,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, true, RVec,
                        true, true,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, false,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        false, true,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, false,
                        true, true, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        false, false, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        false, false, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        false, true, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        false, true, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        true, false, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        true, false, true >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        true, true, false >(
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
         ier = Compute< ClusterIterator, true, MI_OPBC,
                        true, true,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, false,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        false, true,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, false,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, false, Coordinates,
                        true, true,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, false,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        false, true,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, false,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, false, RVec,
                        true, true,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, false,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        false, true,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, false,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, false, MI_OPBC,
                        true, true,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, false,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        false, true,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, false,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, true, Coordinates,
                        true, true,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, false,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        false, true,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, false,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, true, RVec,
                        true, true,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, false,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        false, true,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, false,
                        true, true, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        false, false, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        false, false, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        false, true, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        false, true, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        true, false, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        true, false, true >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        true, true, false >(
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
         ier = Compute< LocatorIterator, true, MI_OPBC,
                        true, true,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, false,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        false, true,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, false,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, false, Coordinates,
                        true, true,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, false,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        false, true,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, false,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, false, RVec,
                        true, true,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, false,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        false, true,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, false,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, false, MI_OPBC,
                        true, true,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, false,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        false, true,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, false,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, true, Coordinates,
                        true, true,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, false,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        false, true,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, false,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, true, RVec,
                        true, true,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, false,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        false, true,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, false,
                        true, true, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        false, false, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        false, false, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        false, true, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        false, true, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        true, false, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        true, false, true >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        true, true, false >(
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
         ier = Compute< IteratorIterator, true, MI_OPBC,
                        true, true,
                        true, true, true >(
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
