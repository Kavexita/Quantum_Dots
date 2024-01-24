!The programm consist in two parts. One create the dataset and the other ones write in a selection

program h5data

  USE hdf5
  USE occup
  IMPLICIT NONE
  
  real(kind=8) :: d(4)
  real(kind=8) :: m(2)
  real(kind=8) :: Ic(30,30)
  real(kind=8) :: emin, emax
  !real(kind=8) :: mu, sigma
  
  CHARACTER(LEN=16), PARAMETER :: filename = "double_dot.h5" ! File name
  CHARACTER(LEN=16), PARAMETER :: mdsetname = "voltages"       ! Dataset name
  CHARACTER(LEN=16), PARAMETER :: vdsetname = "gates_w"
  CHARACTER(LEN=16), PARAMETER :: groupname = "dataset"
  CHARACTER(LEN=16), PARAMETER :: v1dsetname = "V1"
  CHARACTER(LEN=16), PARAMETER :: v2dsetname = "V2"
  
  INTEGER(HSIZE_T), PARAMETER :: dset_size = 6000
  INTEGER(HSIZE_T), PARAMETER :: img_size = 30
  INTEGER(HSIZE_T), PARAMETER :: vector_size = 4
  INTEGER(HSIZE_T), PARAMETER :: one = 1 !For avoid kind problems
  INTEGER(HID_T) :: file_id
  INTEGER(HID_T) :: group_id
  INTEGER(HID_T) :: mdataspace, vdataspace, v1dataspace, v2dataspace
  INTEGER(HID_T) :: mdset_id, vdset_id, v1dset_id, v2dset_id
  INTEGER(HID_T) :: mmemspace, vmemspace
  
  INTEGER(HSIZE_T), DIMENSION(1:3) :: dimsm = (/img_size,img_size,one/)
  INTEGER(HSIZE_T), DIMENSION(1:2) :: dimsv = (/vector_size,one/)
  INTEGER(HSIZE_T), DIMENSION(1) :: dimsv1 = (/img_size/)
  INTEGER(HSIZE_T), DIMENSION(1) :: dimsv2 = (/img_size/)
  REAL, DIMENSION(1:img_size,1:img_size,1) :: mdata                  ! Subset buffer
  REAL, DIMENSION(1:vector_size,1) :: vdata
  REAL, DIMENSION(1:img_size) :: v1data
  REAL, DIMENSION(1:img_size) :: v2data
  INTEGER :: dim0_sub = img_size
  INTEGER :: dim1_sub = img_size
  INTEGER :: dim2_sub = 1
  
  INTEGER(HSIZE_T), DIMENSION(1:3) :: mcount = (/img_size,img_size,one/)  ! Size of hyperslab
  INTEGER(HSIZE_T), DIMENSION(1:3) :: moffset ! Hyperslab offset
  INTEGER(HSIZE_T), DIMENSION(1:3) :: mstride = (/1,1,1/) ! Hyperslab stride
  INTEGER(HSIZE_T), DIMENSION(1:3) :: mblock = (/1,1,1/)  ! Hyperslab block size
  
  INTEGER(HSIZE_T), DIMENSION(1:2) :: vcount = (/vector_size,one/)  ! Size of hyperslab
  INTEGER(HSIZE_T), DIMENSION(1:2) :: voffset ! Hyperslab offset
  INTEGER(HSIZE_T), DIMENSION(1:2) :: vstride = (/1,1/) ! Hyperslab stride
  INTEGER(HSIZE_T), DIMENSION(1:2) :: vblock = (/1,1/)  ! Hyperslab block size
  
  INTEGER(HSIZE_T), DIMENSION(1:3) :: mdimsf = (/img_size,img_size,dset_size/) ! Dataset dimensions
  INTEGER(HSIZE_T), DIMENSION(1:2) :: vdimsf = (/vector_size,dset_size/)
  
  INTEGER :: mrank = 3      ! Dataset rank ( in file )
  INTEGER :: vrank = 2
  INTEGER :: v1rank = 1
  INTEGER :: v2rank = 1
  INTEGER :: error, i, j, k
  
  REAL(HSIZE_T), DIMENSION(1:2) :: mu = (/1.,0./)
  REAL(HSIZE_T), DIMENSION(1:2) :: sigma = (/0.1,0.1/)
  
  !Open the interface
  CALL h5open_f(error)
  
  !Create the file
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
  
  
  !Write the voltajes of the gates
  emin=-40
  emax=+20
  
  do j=1,img_size
  	v1data(j) = emin + ((j-1)/dble(img_size-1))*(emax-emin)
  	v2data(j) = emin + ((j-1)/dble(img_size-1))*(emax-emin)
  enddo
  
  CALL h5screate_simple_f(v1rank, dimsv1, v1dataspace, error)
  
  CALL h5dcreate_f(file_id, v1dsetname, H5T_NATIVE_REAL, v1dataspace, &
       v1dset_id, error)
  
  CALL h5dwrite_f(v1dset_id, H5T_NATIVE_REAL, v1data, dimsv1, error)
  
  CALL h5screate_simple_f(v2rank, dimsv2, v2dataspace, error)
  
  CALL h5dcreate_f(file_id, v2dsetname, H5T_NATIVE_REAL, v2dataspace, &
       v2dset_id, error)
  
  CALL h5dwrite_f(v2dset_id, H5T_NATIVE_REAL, v2data, dimsv2, error)
  
  
  !Create the group
  CALL h5gcreate_f(file_id, groupname, group_id, error)

  !Create the memory space for the matrix dataset
  CALL h5screate_simple_f(mrank, mdimsf, mdataspace, error)
  
  !Create the matrix dataset
  CALL h5dcreate_f(group_id, mdsetname, H5T_NATIVE_REAL, mdataspace, &
       mdset_id, error)
  
  !Create the memory space for the vector dataset
  CALL h5screate_simple_f(vrank, vdimsf, vdataspace, error)
  
  !Create the vector dataset
  CALL h5dcreate_f(group_id, vdsetname, H5T_NATIVE_REAL, vdataspace, &
       vdset_id, error)
  
  !Close the dataspaces, datasets and file
  CALL h5sclose_f(mdataspace, error)
  CALL h5dclose_f(mdset_id, error)
  CALL h5sclose_f(vdataspace, error)
  CALL h5sclose_f(v1dataspace, error)
  CALL h5sclose_f(v2dataspace, error)
  CALL h5dclose_f(vdset_id, error)
  CALL h5dclose_f(v1dset_id, error)
  CALL h5dclose_f(v2dset_id, error)
  CALL h5fclose_f(file_id, error)
  
  CALL h5close_f(error)
  
  !
  !DATASET AND GROUP CREATED AND CLOSED
  !
  
  CALL h5open_f(error)
  
  CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
  
  
  do i = 0,dset_size-1
  
	  call normal_distribution_example(mu(1), sigma(1), d(1))
	  call normal_distribution_example(mu(1), sigma(1), d(4))

	  call normal_distribution_example(mu(2), sigma(2), d(2))
	  call normal_distribution_example(mu(2), sigma(2), d(3))
	  d(2)= abs(d(2))
	  d(3)= abs(d(3))
	  
	  !write(*,*) d

	  ! Call the function
	  call ocupaciones(d, Ic)
	  
	  !Write this values in a file directly in the subroutine
	  
	  moffset = (/0,0,i/)
	  voffset = (/0,i/)
	  
	  mdata(1:img_size,1:img_size,1) = Ic(1:img_size,1:img_size)
          vdata(1:vector_size,1) = d(1:vector_size)
          
          !write(*,*) vdata
	  
	  CALL h5dopen_f(group_id, mdsetname, mdset_id, error)
	  
	  CALL h5dget_space_f(mdset_id, mdataspace, error)
	  CALL h5sselect_hyperslab_f(mdataspace, H5S_SELECT_SET_F, &
	       moffset, mcount, error, mstride, mblock)
	  
	  CALL h5screate_simple_f(mrank, dimsm, mmemspace, error)

	  
	  !CALL h5dwrite_f(mdset_id, H5T_NATIVE_REAL, mdata, dimsm, error, &
	       !mmemspace, mdataspace)
	  
	  CALL h5dwrite_f(mdset_id, H5T_NATIVE_REAL, mdata, dimsm, error, &
	       mmemspace, mdataspace)
	  
	  CALL h5sclose_f(mdataspace, error)
	  CALL h5sclose_f(mmemspace, error)
	  CALL h5dclose_f(mdset_id, error)
	  
	  !Here vector
	  CALL h5dopen_f(group_id, vdsetname, vdset_id, error)
	  
	  CALL h5dget_space_f(vdset_id, vdataspace, error)
	  CALL h5sselect_hyperslab_f(vdataspace, H5S_SELECT_SET_F, &
	       voffset, vcount, error, vstride, vblock)
	  
	  CALL h5screate_simple_f(vrank, dimsv, vmemspace, error)

	  
	  !CALL h5dwrite_f(vdset_id, H5T_NATIVE_REAL, vdata, dimsv, error, &
	       !vmemspace, vdataspace)
	  
	  CALL h5dwrite_f(vdset_id, H5T_NATIVE_REAL, vdata, dimsv, error, &
	       vmemspace, vdataspace)
	  
	  CALL h5sclose_f(vdataspace, error)
	  CALL h5sclose_f(vmemspace, error)
	  CALL h5dclose_f(vdset_id, error)
	  
  enddo
  
  CALL h5fclose_f(file_id, error)
  CALL h5close_f(error)
  
end program h5data
