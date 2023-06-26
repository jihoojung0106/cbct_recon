1. 준비

1) 폴더 구조
-projection data들을 main.c와 같은 경로에 생성된 "input"이라는 folder에 넣습니다.
-같은 폴더에 output을 담을 output 폴더를 만듭니다.

└─ MyProject/
   ├─ input/
   ├─ output/
   ├─ main.c
   └─ Makefile

2) library 설치
-fftw 라이브러리(https://www.fftw.org/) 설치하기


2. 실행

-make run


3. 결과 확인

-imagej 이용해서 750x750, 16bit signed Int로 확인합니다
