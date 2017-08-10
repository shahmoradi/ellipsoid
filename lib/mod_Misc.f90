module Misc
  
  implicit none
  
  integer, parameter, private :: maxFileRecordLength = 2047
  
  type DynamicString
    character (len=:), allocatable :: record
  end type DynamicString

  interface num2str
    module procedure int2str
  end interface

contains

!*********************************************************************
!*********************************************************************

  recursive function replaceStr(string,search,substitute) result(modifiedString)
    implicit none
    character(len=*), intent(in)  :: string, search, substitute
    character(len=:), allocatable :: modifiedString
    integer                       :: i, stringLen, searchLen,currentPos
    stringLen = len(string)
    searchLen = len(search)
    if (stringLen<searchLen) then
      modifiedString = string
      return
    end if
    i = 1
    do
      if (string(i:i+searchLen-1)==search) then
        modifiedString = string(1:i-1) // substitute // replaceStr(string(i+searchLen:stringLen),search,substitute)
        exit
      end if
      if (i+searchLen>stringLen) then
        modifiedString = string
        exit
      end if
      i = i + 1
      cycle
    end do
  end function replaceStr

!*********************************************************************
!*********************************************************************

  function makeFormat(count,type,len,precision)
    implicit none
    integer     , intent(in)  :: count,len,precision
    character(1), intent(in)  :: type
    character(:), allocatable :: makeFormat
    makeFormat = '('//num2str(count)//type//num2str(len)//'.'//num2str(precision)//')'
  end function makeFormat

!*********************************************************************
!*********************************************************************

  pure function getLowerCase(string)
    ! Changes a string to upper case
    implicit None
    character(*), intent(in) :: string
    character(len(string))   :: getLowerCase
    character(26), parameter :: lowerCase = 'abcdefghijklmnopqrstuvwxyz', upperCase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer                  :: ic, i
    ! capitalize each letter if it is lowercase
    getLowerCase = string
    do i = 1, len(string)
        ic = INDEX(upperCase, string(i:i))
        if (ic > 0) getLowerCase(i:i) = lowerCase(ic:ic)
    end do
  end function getLowerCase

!*********************************************************************
!*********************************************************************

  pure function getUpperCase(string)
    ! Changes a string to upper case
    implicit None
    character(*), intent(in) :: string
    character(len(string))   :: getUpperCase
    character(26), parameter :: lowerCase = 'abcdefghijklmnopqrstuvwxyz', upperCase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer                  :: ic, i
    ! capitalize each letter if it is lowercase
    getUpperCase = string
    do i = 1, len(string)
        ic = INDEX(lowerCase, string(i:i))
        if (ic > 0) getUpperCase(i:i) = upperCase(ic:ic)
    end do
  end function getUpperCase

!*********************************************************************
!*********************************************************************

  function getFileExt(fileName)
    implicit none
    character(len=*), intent(in)  :: fileName
    character(len=:), allocatable :: getFileExt
    integer                          :: dotPos, lenFileName
    lenFileName = len(fileName)
    dotPos = scan(fileName,'.',back=.true.)
    if (dotPos>0 .and. dotPos<lenFileName) then
      getFileExt = fileName(dotPos+1:lenFileName)
    else
      getFileExt = ""
    end if
  end function getFileExt

!*********************************************************************
!*********************************************************************

  subroutine getFileBaseExt(fileName,fileBase,fileExt)
    implicit none
    character(len=*)             , intent(in)  :: fileName
    character(len=:), allocatable, intent(out) :: fileBase, fileExt
    integer                          :: dotPos, lenFileName
    lenFileName = len(fileName)
    dotPos = scan(fileName,'.',back=.true.)
    if (dotPos<1 .or. dotPos>lenFileName) then
      fileBase = fileName
      fileExt  = ""
    elseif (dotPos==1) then
      fileBase = ""
      fileExt  = fileName(2:lenFileName)
    elseif (dotPos==lenFileName) then
      fileBase = fileName(1:lenFileName-1)
      fileExt  = ""
    else
      fileBase = fileName(1:dotPos-1)
      fileExt  = fileName(dotPos+1:lenFileName)
    end if
  end subroutine getFileBaseExt

!*********************************************************************
!*********************************************************************

  !function SplitStr(string,separationChar)
  !  implicit none
  !  character(*), intent(in)         :: string, separationChar
  !  type(DynamicString), allocatable :: SplitStr, SplitStrDummy
  !  integer                          :: stringLen, i, separationLen
  !  character(2047)                  :: record
  !  stringLen = len(string)
  !  separationLen = len(string)
  !  if (separationChar>=string) return
  !  i = 1
  !  do
  !    
  !    if (i==stringLen-separationLen) exit
  !    i = i + 1
  !    cycle
  !  end do
  !end function SplitStr

!*********************************************************************
!*********************************************************************

  function getOS()
    implicit none
    character(len=:)  , allocatable :: getOS
    character(len=:)  , allocatable :: filename
    character(len=7)                :: os
    integer                         :: idummy
    character(8)                    :: date
    character(10)                   :: time
    character(10)                   :: record
    logical                         :: exist
    call get_environment_variable('OS',os)
    if (os=='Windows') then
      getOS = 'Windows'
    else    ! it's either Linux- or Darwin- based OS
      idummy = 0
      do
        idummy = idummy + 1
        call date_and_time(date,time)
        filename = date // '_' // time // '_' // 'getOS_' // num2str(idummy) // '.temp'
        inquire(file=filename,exist=exist)    ! check if the file already exists
        if (exist) cycle
        exit
      end do
      call execute_command_line('uname > '//filename)
      open(newunit=idummy,file=filename,status='old')
      read(idummy,*) os
      close(idummy)
      call execute_command_line('rm '//filename)
      getOS = trim(adjustl(os))
    end if
  end function getOS

!*********************************************************************
!*********************************************************************

  function getFileList(searchText,order,excludeText)

    implicit none
    character(len=*)   , intent(in)           :: searchText
    character(len=*)   , intent(in), optional :: order
    character(len=*)   , intent(in), optional :: excludeText
    type(DynamicString), allocatable          :: getFileList(:)
    character(len=:)   , allocatable          :: command,filename,orderLowerCase
    character(len=maxFileRecordLength)        :: record
    integer                                   :: iunit,counter,iostat,nRecord,nskip
    character(8)                              :: date
    character(10)                             :: time
    logical                                   :: exist
    
    if (present(order)) then
      orderLowerCase = getLowerCase(order)
    else
      orderLowerCase = 'name'
    end if

    if (getSlash()=='\') then  ! it's Windows cmd
      if (orderLowerCase=='name') then  ! ascending in name
        command = 'dir /b /a-d ' // searchText
      elseif (orderLowerCase=='date') then   ! oldest will be first
        command = 'dir /b /a-d /o:d ' // searchText
      else
        write(*,*) '    FATAL: In Misc@getFileList()'
        write(*,*) '           The requested file order is not supported.'
        write(*,*) '           order = ', orderLowerCase
        write(*,*) 'Program aborted.'
      end if
      if (present(excludeText)) then
        command = command // " | findstr /v /i " // trim(adjustl(excludeText))
      end if
    else
      if (orderLowerCase=='name') then  ! ascending in name
        command = 'ls -1 ' // searchText
      elseif (orderLowerCase=='date') then   ! oldest will be first
        command = 'ls -tr ' // searchText
      else
        write(*,*) '    FATAL: In Misc@getFileList()'
        write(*,*) '           The requested file order is not supported.'
        write(*,*) '           order = ', orderLowerCase
        write(*,*) 'Program aborted.'
      end if
      if (present(excludeText)) then
        command = command // " --ignore=" // trim(adjustl(excludeText))
      end if
    end if

    ! generate a brand new, non-existing filename
    counter = 0
    do
      counter = counter + 1
      call date_and_time(date,time)
      filename = date // '_' // time // '_' // 'getFileList_' // num2str(counter) // '.temp'
      inquire(file=filename,exist=exist)    ! check if the file already exists
      if (exist) cycle
      exit
    end do
    call execute_command_line(command//' > '//filename)

    nRecord = getNumRecordInFile(filename)
    
    ! check filename is not among records
    nskip = 0
    open(newunit=iunit,file=filename,status='old')
    do counter = 1,nRecord
      read(iunit,'(A)',iostat=iostat) record
      if(iostat==0) then
        if(filename==trim(adjustl(record))) nskip = nskip + 1
      else
        write(*,*) '    FATAL (1): In Misc@getFileList()'
        write(*,*) '               Error occurred while reading file.'
        write(*,*) 'Program aborted.'
        stop
      end if
    end do
    close(iunit)

    allocate(getFileList(nRecord-nskip))
    open(newunit=iunit,file=filename,status='old')
    do counter = 1,nRecord
      read(iunit,'(A)',iostat=iostat) record
      if(iostat==0) then
        if (filename/=trim(adjustl(record))) getFileList(counter)%record = trim(adjustl(record))
      else
        write(*,*) '    FATAL (2): In Misc@getFileList()'
        write(*,*) '               Error occurred while reading file.'
        write(*,*) 'Program aborted.'
        stop
      end if
    end do
    close(iunit)
    
    if (getSlash()=='\') then  ! it's Windows cmd
      command = 'del '//filename
    else
      command = 'rm '//filename
    end if
    call execute_command_line(command)

  end function getFileList

!*********************************************************************
!*********************************************************************

  function getSystemInfo()

    implicit none
    type(DynamicString), allocatable   :: getSystemInfo(:)
    character(len=:)   , allocatable   :: command,filename
    character(len=maxFileRecordLength) :: record
    integer                            :: iunit,counter,iostat,nRecord
    character(8)                       :: date
    character(10)                      :: time
    logical                            :: exist

    ! generate a brand new, non-existing filename
    counter = 0
    do
      counter = counter + 1
      call date_and_time(date,time)
      filename = date // '_' // time // '_' // 'getSystemInfo_' // num2str(counter) // '.temp'
      inquire(file=filename,exist=exist)    ! check if the file already exists
      if (exist) cycle
      exit
    end do    
    
    if (getSlash()=='\') then  ! it's Windows cmd
      command = 'systeminfo > '//filename
    else
      command = 'uname -a > '//filename//'; sudo lshw -short > '//filename//'; lscpu > '//filename
    end if

    call execute_command_line(command)

    nRecord = getNumRecordInFile(filename)
    
    allocate(getSystemInfo(nRecord))
    open(newunit=iunit,file=filename,status='old')
    do counter = 1,nRecord
      read(iunit,'(A)',iostat=iostat) record
      if(iostat==0) then
        getSystemInfo(counter)%record = trim(adjustl(record))
      else
        write(*,*) 'FATAL error occurred while reading file in Misc in getSystemInfo().'
        write(*,*) 'Program aborted.'
        stop
      end if
    end do
    close(iunit)
    
    if (getSlash()=='\') then  ! it's Windows cmd
      command = 'del '//filename
    else
      command = 'rm '//filename
    end if
    call execute_command_line(command)

  end function getSystemInfo

!*********************************************************************
!*********************************************************************

  function getFileContent(filename)
    implicit none
    character(len=*)    , intent(in)   :: filename
    type(DynamicString) , allocatable  :: getFileContent(:)
    character(len=maxFileRecordLength) :: record
    integer                            :: nRecord,iunit,irecord,iostat
    nRecord = getNumRecordInFile(filename)
    allocate(getFileContent(nRecord))
    open(newunit=iunit,file=filename,status='old')
    do irecord = 1,nRecord
      read(iunit,'(A)',iostat=iostat) record
      if(iostat==0) then
        getFileContent(irecord)%record = trim(adjustl(record))
      else
        write(*,*) 'FATAL error occurred while reading file in Misc in getFileContent().'
        write(*,*) 'Program aborted.'
        stop
      end if
    end do
    close(iunit)
  end function getFileContent

!*********************************************************************
!*********************************************************************

  function getNumRecordInFile(filename)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=8)             :: record
    integer                      :: getNumRecordInFile,iunit,iostat
    open(newunit=iunit,file=filename,status='old')
    getNumRecordInFile = 0
    do
      read(iunit,'(A)',iostat=iostat) record
      if(iostat==0) then
        getNumRecordInFile = getNumRecordInFile + 1
        cycle
      elseif(iostat<0) then
        exit
      else
        write(*,*) 'FATAL error occurred reading file in Misc in getNumRecordInFile().'
        write(*,*) 'Program aborted.'
        stop
      end if
    end do
    close(iunit)
  end function getNumRecordInFile

!*********************************************************************
!*********************************************************************

  subroutine splitPath(path,directory,filename)
    implicit none
    character(len=*)             , intent(in)  :: path
    character(len=:), allocatable, intent(out) :: directory, filename
    character(len=:), allocatable              :: slash,pathDummy
    integer                                    :: pathLength, slashPosition
    slash = getSlash()
    if (slash=='\') then
      pathDummy = replaceStr(trim(adjustl(path)),'/','\')
    else
      pathDummy = trim(adjustl(path))
    end if
    pathLength = len(pathDummy)
    slashPosition = pathLength
    do
      if (pathDummy(slashPosition:slashPosition)==slash .or. slashPosition==1) exit
      slashPosition = slashPosition - 1
      cycle
    end do
    if (slashPosition==pathLength) then
      filename  = ''
      directory = pathDummy
    elseif (slashPosition==1) then
      filename  = pathDummy
      directory = ''
    else
      filename  = pathDummy(slashPosition+1:pathLength)
      directory = pathDummy(1:slashPosition)
    end if
  end subroutine splitPath

!*********************************************************************
!*********************************************************************

  character(len=1) function getSlash()
    implicit none
    character(len=7) :: os
    call get_environment_variable('OS',os)
    if (os=='Windows') then
      getSlash = '\'
    else
      getSlash = '/'
    end if
  end function getSlash

!*********************************************************************
!*********************************************************************
   
  function getNiceDateTime()
    implicit none
    character(len=21)             :: getNiceDateTime
    character(10)                 :: time
    character(8)                  :: date
    call date_and_time(date,time)
    getNiceDateTime = date(1:4)//'/'//date(5:6)//'/'//date(7:8)//' - '//time(1:2)//':'//time(3:4)//':'//time(5:6)
  end function getNiceDateTime

!*********************************************************************
!*********************************************************************
   
  pure function int2str(integerIn,formatIn)
    implicit none
    integer     , intent(in)           :: integerIn
    character(*), intent(in), optional :: formatIn
    character(:), allocatable          :: int2str
    integer                            :: i,length
    character(len=63)                  :: thisFormat
    if (present(formatIn)) then
      write(thisFormat,formatIn) integerIn
      int2str = trim(adjustl(thisFormat))
    else
      do i=1,63
        if(abs(integerIn)<10**i) then
          length = i
          if (integerIn<0) length = length + 1
          exit
        end if
      end do
      allocate(character(length) :: int2str)
      write(thisFormat,'(1I63)') length
      thisFormat = '(1I' // trim(adjustl(thisFormat)) // ')'
      write(int2str,thisFormat) integerIn
    end if
  end function int2str

!*********************************************************************
!*********************************************************************

  integer elemental function str2int(strIn)
    implicit none
    character(len=*), intent(in) :: strIn
    integer                      :: iostat
    character(len=63)            :: thisFormat
    read(strIn,*,iostat=iostat) str2int
    !if (iostat/=0) then
    !  write(*,*) 'FATAL: str2int() in module Misc failed.'
    !  stop
    !end if
  end function str2int

!*********************************************************************
!*********************************************************************  
  
end module Misc