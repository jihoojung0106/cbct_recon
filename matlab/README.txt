1. 준비

1)airmap

-airmap 파일을 airmap_fit.raw 로 이름 지정하여 main.m과 같은 경로에 넣습니다.

2) 폴더 구조
-projection data들을 main.m와 같은 경로에 생성된 "input"이라는 folder에 넣습니다.
-같은 폴더에 output을 담을 output 폴더를 만듭니다.

└─ MyProject/
   ├─ input/
   ├─ output/
   ├─ airmap_fit.raw
   ├─ main.m
   ├─ MyParameter.m
   ├─ filtering.m
   ├─ CTbackprojection.m
   ├─ inpaint_nans.m		
   └─ backprojection.m


3) toolbox

- Parallel Computing Toolbox(https://www.mathworks.com/products/parallel-computing.html) 설치합니다.


2. 실행

-매트랩에서 실행


3. 결과 확인

-image j, 750x750 16 bit int 로 확인합니다