module types
use const
implicit none


type :: t_block
      integer(INT_KIND) :: nCells(3)
      integer(INT_KIND) :: nPoints(3)
      real(REAL_KIND), allocatable :: coords(:,:,:,:)
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


      integer(INT_KIND), allocatable :: edge_points(:,:)
      !< Indices of points at edge end (point_id, edge_id)

      real(REAL_KIND), allocatable   :: edge_lengths(:)

      real(REAL_KIND), allocatable   :: edge_springs(:)

      real(REAL_KIND), allocatable   :: edge_vectors(:,:)

      real(REAL_KIND), allocatable   :: edge_forces(:,:)
      !< FORCE VECTOR of EDGE (edge_spring * edge_vector)

end type t_unstr

end module types
