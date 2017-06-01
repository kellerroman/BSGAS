module types
use const
implicit none


type :: t_boundary_condition
      integer(INT_KIND)              :: bc_type
      !< Boundary-Type: >0: Blocknumber of other Block
      !<                 0: Uninitialized
      integer(INT_KIND)              :: neighbor_face
      integer(INT_KIND)              :: permutation
end type t_boundary_condition

type :: t_block
      integer(INT_KIND)              :: nCells(3)
      integer(INT_KIND)              :: nPoints(3)
      real(REAL_KIND)  , allocatable :: coords(:,:,:,:)
      integer(INT_KIND), allocatable :: refs(:,:,:)
      type(t_boundary_condition), allocatable :: boundary_cond(:)
      !< Boundary Condition at the Block Faces (Face-id)
      !< Face: 1 = West (i = 1), 2 = East (i = Imax)
      !<       3 = South(j = 1), 4 = North(j = Jmax)
      !<       5 = Front(k = 1), 4 = Back (k = kmax)
end type t_block

type :: t_unstr
      integer(INT_KIND) :: dimension
      integer(INT_KIND) :: npoint
      integer(INT_KIND) :: nedge

      real(REAL_KIND), allocatable   :: point_coords(:,:)

      integer(INT_KIND), allocatable :: point_refs(:,:)
      !< Reference Cell in Block struct grid ([b,i,j,k],point_id)

      integer(INT_KIND), allocatable :: point_nedges(:)
      !< Number of edges attached to a point

      integer(INT_KIND), allocatable :: point_edges(:,:)
      !< Indices of attached edges (edge_id,point_id) 

      real(REAL_KIND), allocatable   :: point_edge_signs(:,:)
      !< Edge vectors are always from first to second point.
      !< So for second points, the sign of the force must be inverted

      real(REAL_KIND), allocatable   :: point_forces(:,:)
      !< Force, sum of attateced edge forces, acting on the point

      logical, allocatable           :: point_move_rest(:)
      !< Logical if the Point has Movement restrictions

      real(REAL_KIND), allocatable   :: point_move_rest_vector(:,:)
      !< Vector with the movement restriction


      integer(INT_KIND), allocatable :: edge_points(:,:)
      !< Indices of points at edge end (point_id, edge_id)

      real(REAL_KIND), allocatable   :: edge_lengths(:)

      real(REAL_KIND), allocatable   :: edge_springs(:)

      real(REAL_KIND), allocatable   :: edge_vectors(:,:)

      real(REAL_KIND), allocatable   :: edge_forces(:,:)
      !< FORCE VECTOR of EDGE (edge_spring * edge_vector)

      integer(INT_KIND), allocatable :: edge_nneighbor(:)
      !< NUmber of Neighbor Edges in same direction 
      integer(INT_KIND), allocatable :: edge_neighbor(:,:)
      !< Neighbor Edges in same direction (length diff resitrictions)


      !!!!!! LISTS FOR SPECIAL TYPES:

      integer(INT_KIND) :: nWallEdge
      integer(INT_KIND), allocatable :: wall_edges(:)
      !< List of all edges that are ajectent to a wall thus having size requirements
end type t_unstr

end module types
