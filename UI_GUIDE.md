# ONT Plasmid Analyzer 사용 가이드

이 UI는 외부 서버로 데이터를 전송하지 않고 **현재 Linux 워크스테이션 안에서만**
실행됩니다. 기존 `run_pipeline.sh` 명령줄 방식도 그대로 사용할 수 있습니다.

공개 저장소에는 특정 회사의 이미지 로고를 포함하지 않습니다. 왼쪽 하단에 배경 없는
text wordmark를 사용하려면 Linux server에서 다음을 한 번 실행합니다. 이 설정 파일은
Git에서 자동 제외되며 NASCA로 암호화된 이미지가 필요하지 않습니다.

```bash
./configure_local_branding.sh "SAMSUNG BIOLOGICS" "Jongin Baek"
```

## 1. 설치

기존 `NGS_ONT_env`를 업데이트합니다.

```bash
cd NGS_ONT
conda env update -p /home/MCET03/conda_envs/NGS_ONT_env -f environment.yml
conda activate /home/MCET03/conda_envs/NGS_ONT_env
```

인터넷이 차단된 사내 서버에서는 로컬 미러만 지정할 수 있습니다.

```bash
conda env update -p /home/MCET03/conda_envs/NGS_ONT_env --offline --override-channels \
  -c file:///data/conda_repo/mirror/conda-forge/conda-forge \
  -c file:///data/conda_repo/mirror/bioconda \
  -f environment.yml
```

설치 확인:

```bash
streamlit --version
minimap2 --version
samtools --version
bcftools --version
NanoFilt --version
```

## 2. UI 실행

### 아이콘으로 실행

Linux desktop에서는 최초 한 번만 다음을 실행합니다.

```bash
cd /data/user/MCET03/04_ONT/NGS_ONT
chmod +x launch_ui.sh install_linux_launcher.sh run_ui.sh
./install_linux_launcher.sh
```

이후 바탕화면의 **ONT Plasmid Analyzer** 아이콘을 클릭하면 `NGS_ONT_env` 탐색,
UI 실행 및 브라우저 열기가 자동으로 진행됩니다. 기본 포트는 `8502`입니다.

Windows + WSL 환경에서는 프로젝트의 `Launch_ONT_UI.cmd`를 바탕화면으로 복사한 뒤
더블클릭합니다. 기본 WSL 프로젝트 경로는 `~/NGS_ONT_batch`입니다.

Windows에서 회사 Linux server에 접속할 때는 server에서 최초 한 번 다음을 실행합니다.

```bash
chmod +x ont_one_click.sh launch_ui.sh run_ui.sh run_pipeline.sh
```

그다음 `ONT_Plasmid_Analyzer_One_Click.bat`을 Windows 바탕화면으로 내려받습니다. 이 파일은
MobaXterm 설치 위치를 자동으로 찾고, 기존 tunnel이 정상이면 브라우저만 즉시 엽니다.
Tunnel이 없으면 기존 MobaXterm 창에 SSH tab을 추가하거나 MobaXterm을 새로 시작한 뒤,
server의 `ont_one_click.sh`로 환경 활성화와 Streamlit 재시작을 수행합니다. Windows에서
표시되는 주소는 `http://127.0.0.1:18502`입니다.
열린 SSH 창은 분석 중 닫지 않아야 합니다.

### 명령어로 실행

```bash
conda activate /home/MCET03/conda_envs/NGS_ONT_env
cd /data/user/MCET03/04_ONT/NGS_ONT
./run_ui.sh
```

코드나 branding 설정을 바꾼 뒤 기존 UI를 한 번에 재시작하려면 다음을 사용합니다.

```bash
./launch_ui.sh --restart
```

브라우저에서 다음 주소를 엽니다.

```text
http://<Linux-server-IP>:8501
```

포트를 바꾸려면:

```bash
ONT_UI_PORT=8502 ./run_ui.sh
```

서버 외부 접속을 막고 SSH 포트 포워딩만 사용하려면:

```bash
ONT_UI_ADDRESS=127.0.0.1 ./run_ui.sh
```

## 3. Batch analysis (권장)

먼저 화면 상단에 이번 분석의 `Project name`과 `User name`을 입력합니다. 두 항목은
분석 실행에 필수이며 결과 화면과 server 사용 기록에 남습니다.

1. `Reference 폴더 탐색`에서 폴더 버튼을 눌러 이동하고, FASTA가 들어 있는 폴더에서
   `이 폴더 사용`을 누릅니다. 탐색 시작 위치와 최상위 범위는 `/data`입니다.
2. Reference는 파일명 기준 자연 정렬되며(예: `01, 02, ..., 10`),
   필요할 때만 접힌 `Reference 확인`에서 파일명과 길이를 확인합니다.
3. `ONT 결과 폴더 탐색`에서 sample 하위 폴더가 들어 있는 상위 폴더를 선택합니다.
   표준 ONT run 폴더의 `fastq_pass`만 자동 검색하며 `fastq_fail`, `other_reports`,
   `pod5_*`는 제외합니다. FASTQ/FASTQ.GZ 원본은 복사하지 않고 직접 연결합니다.
4. `Reference별 ONT sample 선택`에서 분석할 조합을 직접 지정합니다. 이름이나 순서를
   기준으로 자동 배정하지 않습니다. Reference마다 sample 수를 다르게 선택할 수 있고,
   같은 sample을 여러 Reference와 비교할 수도 있습니다. 선택하지 않은 Reference와
   ONT sample은 이번 분석에서 제외됩니다.
5. 왼쪽 하단의 `Analysis settings` 버튼을 누르면 중앙에 나타나는 설정 화면에서
   CPU, 병렬 sample 수, read filter 및 variant review 기준을 저장할 수 있습니다.
6. 중복·누락 경고가 없는지 확인한 뒤 `Batch analysis 실행`을 누릅니다.

`분석 정보`, `Reference 입력`, `ONT sample 입력`, `Reference별 ONT sample 선택`은
각각 별도의 border 영역으로 구분되어 있습니다.

왼쪽 sidebar의 작은 `분석 도구` checklist에서 `minimap2`, `samtools`, `bcftools`,
`gzip`, `bash`가 모두 ✅인지 확인할 수 있습니다. Batch QC는 설치된 `chopper`를 우선 사용하며,
없으면 `NanoFilt`를 사용합니다. 둘 다 없으면 로그에 필터 생략을 표시합니다.
FASTQ 압축에는 `pigz`가 있으면 이를 사용하고, 없으면 `gzip`을 사용합니다. Sidebar를 접어도 화면 왼쪽 위의 화살표로 다시
열 수 있습니다.

각 입력 항목의 `?`에 커서를 올리면 설정 의미와 권장 사용법을 확인할 수 있습니다.

Primary bcftools VCF에 call이 없지만 samtools consensus와 Reference 사이에 차이가
있으면 mutation을 숨기지 않고 `REVIEW`로 표시합니다. 이 경우 상세 결과에
consensus 기반 결과라는 안내가 함께 표시되며, 실제 판정 전에는 원본 read와 depth를
확인하는 것이 좋습니다.

분석 완료 후 `Intermediate QC 및 결과 파일`을 열고 sample 하나를 선택하면 1/6~6/6
완료 상태뿐 아니라 단계별 핵심 결과도 확인할 수 있습니다. Merge/QC 단계에는 read 수,
Read N50, mean read Q, Q≥30 read 비율과 filter 통과율이 표시되고, Mapping 단계에는
mapping rate, Mean depth, ≥10× coverage가 표시됩니다. Raw consensus는 길이와 N base 수,
Variant calling은 SNP/INS/DEL 및 PASS/REVIEW 수, Final result는 consensus 길이와 최종
판정을 보여줍니다. 이 값은 `<sample>.fastq_qc.json`에 한 번 계산해 저장하므로 UI를
새로고침할 때 대용량 FASTQ를 다시 읽지 않습니다.
Raw/Final consensus, Primary VCF, QC report는 바로 받을 수 있습니다. FASTQ, BAM/index,
depth 등 전체 중간 산출물은 `Complete results ZIP 준비`를 눌렀을 때만 묶이므로 평소에는
대규모 batch의 UI 속도에 영향을 주지 않습니다.

Offline 환경에서 `pyarrow`와 `arrow-cpp`가 호환되지 않으면 UI는 자동으로 Markdown
table을 사용합니다. Pyarrow가 정상 import되는 환경에서는 interactive table과 depth
chart가 자동 활성화되며, 어느 경우에도 variant calling 결과에는 차이가 없습니다.

결과는 reference별로 연결된 sample이 묶여 표시되며 `CLEAN`,
`VARIANT DETECTED`, `REVIEW`, `ERROR` 상태를 제공합니다. 전체 결과 요약은 CSV로
내려받을 수 있고 BAM, VCF, consensus FASTA 및 보고서는 서버 결과 폴더에 남습니다.
실행 command, 상세 log, 서버 결과 경로는 UI에 표시되지 않습니다. 관리자 진단을 위한
log 파일은 서버 실행 폴더에 계속 보존됩니다.

별도의 command 출력 없이 다음 사용 기록이 server에 저장됩니다.

- `usage_logs/analysis_usage.csv`: KST timestamp, STARTED/COMPLETED/FAILED,
  Project name, User name, Reference 수, 분석 조합 수, job ID
- `usage_logs/ui_access.log`: `.bat`을 통해 UI를 실제 시작한 timestamp와 Windows 사용자명

Sample ID는 `barcode13` 같은 번호, ONT sample alias 폴더명 또는 FASTQ 파일명에서
자동으로 가져옵니다. 따라서 barcode라는 이름을 반드시 사용할 필요가 없습니다.

> 기본 탐색 범위는 `/data`입니다. 다른 범위가 꼭 필요한 경우에만 UI 실행 전에
> `ONT_SERVER_ROOT` 환경변수를 지정합니다.

## 4. 결과 위치

각 실행은 서로 섞이지 않게 다음 경로에 저장됩니다.

```text
ui_runs/<실행시각>_<experiment>_<job-id>/
├── config.yaml
├── pipeline.log
├── references/
├── data/          # 서버 FASTQ는 복사하지 않고 symbolic link로 연결
└── results/
```

UI에서 CSV, VCF.GZ, consensus FASTA, sample report를 내려받을 수 있습니다.
BAM/BAM.BAI와 depth 파일은 크기가 클 수 있어 결과 폴더에만 보관합니다.

## 5. 정확도 해석

- SNP는 이 UI에서 1-bp substitution/point mutation을 의미합니다.
- `PASS`는 현재 선택한 QUAL, DP, allele-fraction, edge 기준을 만족합니다.
- `REVIEW`는 오류 확정이 아니라 추가 확인이 필요한 call입니다.
- ONT indel은 homopolymer에서 false positive가 증가할 수 있으므로 BAM을 IGV로
  열어 read-level support를 함께 확인하는 것이 좋습니다.
- plasmid/vector 단일 clone의 진짜 변이는 보통 높은 allele fraction을 보입니다.

## 6. 성능 최적화

- Batch pipeline을 시작하면 서버가 각 실행 폴더의 `pipeline.log`에 직접 기록합니다. 브라우저 연결이 끊겨도 시작된 pipeline은 계속 실행됩니다.
- QC 통계는 전체 입력과 필터 통과 read에서 각각 계산합니다. Screening depth는 매핑에 사용할 read 선택에 적용됩니다.
- 공유 filesystem에서 불필요한 source scan이 발생하지 않도록 Streamlit file watcher를 끕니다.
- 진행 화면은 완료된 sample 수만 갱신합니다.
- Sample ID 편집표와 reference별 상세 결과는 사용자가 요청할 때만 생성합니다.
- Analysis log는 사용자가 선택할 때 최근 300줄만 표시합니다.

## 7. 테스트

```bash
python -m unittest discover -s tests -v
bash -n run_pipeline.sh run_ui.sh launch_ui.sh install_linux_launcher.sh
```

## 실행 중 분석 재연결 · Read QC 진단

- 분석은 서버에서 브라우저와 분리되어 실행됩니다. 창을 닫거나 새로고침해도 계속되며,
  다시 열면 "서버에서 분석이 실행 중입니다" 안내와 로그 끝부분, **완료까지 따라가기 / 이 분석 중지**
  버튼이 표시됩니다. 실행 중에는 새 분석 실행 버튼이 잠겨 같은 분석이 중복 실행되지 않습니다.
- 이전 세션에서 끝난 분석은 **최근 완료된 분석 결과 불러오기**로 다시 열 수 있습니다.
- `threads × parallel_jobs`가 CPU 코어 수를 넘으면 자동으로 줄여서 실행합니다.
- Read QC에 **Basecall model**(FASTQ header의 Dorado 모델), **Full-length reads**(reference 길이 90-110% read
  비율), **N50 / reference**, read 길이 분포 그래프, 그리고 해석 힌트가 표시됩니다.
  - HAC/FAST 모델이면 Mean read Q는 보통 Q15-18입니다. Q20 이상이 필요하면 POD5를 SUP 모델로 다시 basecalling하세요.
  - N50이 reference의 80% 미만이면 DNA 절단·분해, 복수 제한효소 자리, 짧은 오염 DNA를 의심합니다.
- Read 필터(길이/mean Q)는 `ont_ui/fastq_qc.py`가 QC와 같은 한 번의 읽기에서 수행합니다
  (NanoFilt/chopper와 같은 mean-Q 정의, NanoFilt 결과와 동일한 출력을 확인함). NanoFilt 설치는 더 이상 필요 없습니다.

## CPU 자동 배분

`Analysis settings`에서 Threads / Parallel samples가 0(기본값)이면 자동 배분합니다.
샘플이 많으면 샘플당 2 threads로 더 많은 샘플을 동시에 처리하고(256코어·96샘플 → 2 × 96),
샘플이 적으면 샘플당 최대 16 threads를 씁니다(256코어·5샘플 → 16 × 5).

## Dorado basecalling (GPU)

사이드바 **⚡ Dorado basecalling (GPU)** 화면에서 POD5를 다시 basecalling합니다.

1. 인터넷 PC에서 모델 다운로드 후 서버의 모델 폴더(기본 `/data/dorado_models`, `ONT_DORADO_MODELS`로 변경)에 복사:
   `dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.2.0 --models-directory <폴더>`
2. POD5 폴더 선택 → 모델(SUP 권장) → Kit(SQK-RBK114-96 / SQK-NBD114-96 등) → **Basecalling 시작**.
3. 결과는 `<출력 폴더>/fastq_pass/barcodeNN/*.fastq.gz`로 정리됩니다. Batch analysis에서 이 `fastq_pass` 폴더를 FASTQ로 선택하세요.
4. 작업은 브라우저와 분리되어 실행되며, 같은 화면에서 진행 로그 확인·중지가 가능합니다. GPU 작업은 한 번에 하나만 실행합니다.

Dorado 실행 파일은 PATH 또는 `ONT_DORADO`로 지정합니다. Blackwell GPU는 이를 지원하는 최신 Dorado 버전이 필요합니다.

## medaka variant caller

`Analysis settings` → **Variant caller: medaka**. ONT 신경망 기반 호출로 homopolymer 주변 indel 판정이 bcftools보다 정확합니다.
- 모델을 비워두면 FASTQ header의 basecall 모델로 자동 선택합니다(HAC v5.x → `r1041_e82_400bps_hac_variant_v5.0.0`).
  pip 설치 medaka 2.x에는 r10.4.1 400bps HAC/SUP variant 모델이 포함되어 오프라인에서도 동작합니다.
- `medaka device: gpu`면 CUDA GPU를 사용하고, `medaka GPU 동시 실행 수`로 동시에 GPU를 쓰는 샘플 수를 제한합니다.
- medaka가 없거나 실패하면 자동으로 bcftools로 대체하고 로그에 WARNING을 남깁니다.
- medaka 2.x(`medaka_variant`)와 1.x(`medaka_haploid_variant`)를 모두 지원합니다. `bgzip`/`tabix`(htslib)가 필요합니다.
