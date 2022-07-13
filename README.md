# EP3_Numerico

# Numerico_EP1

## Requisitos
* Python 3.x
* Numpy
* Math

## Conteúdo
O arquivo funcoes.py contém as funções que foram utilizadas para resolver o problema
proposto. Cada função tem seu comportamento explicado no código. O arquivo main.py é o 
arquivo que deve ser executado para obter os resultados do programa. Informações mais
detalhadas estão na parte execução desse LEIAME.

| Entrada  |   Tipo    | 
|----------|-----------|
| n        |   `int`   |
| k(x)     |   `str`   | 
| q(x)     |   `str`   | 
| L        |   `int`   | 
| u(0)     |  `float`  | 
| u(L)     |  `float`  | 

| Saída    | Tipo       | 
|----------|------------|
| Vetor x  | `np.array` | 
| Vetor y  | `np.array` | 



## Execução
### Iniciando o programa
Para executar o programa, você possui duas opções. Pode utilizar a função RUN de um 
interpretador de sua escolha ou pode rodar através do terminal, utilizando python 3.x. Para 
tal, basta que você primeiro navegue até o local do arquivo
`cd local\do\arquivo`
E depois execute-o digitando
`python3 main.py`
Feito isso, o programa terá sido inciado e você deverá está vendo uma mensagem de boas vindas.

### Modo de operação
Antes de inserir dados você devera especificar qual o modo de operação desejado e tens três 
escolhas: 
    [1] - MÉTODO DE RITZ RALEIGH CÚBICO
    [2] - MÉTODO DAS DIFERENÇAS FINITAS LINEAR
    [3] - MÉTODO DOS ELEMENTOS FINITOS
Para selecionar uma opção digite APENAS o número em colchetes e pressione enter.


### Modo de inserção de dados
Após especificar o modo de operação serão pedidos os dados relevantes para a resolução
do problema de maneira sequencial. 

    INSERÇÃO DE FUNÇÕES

    Atente-se ao modo de inserir equações. Sempre insira a equação entre parentesis e 
    seguindo as regras abaixo:

    | OPERAÇÃO             | OPERADOR  | 
    |----------------------|-----------|
    | Multiplicação        |     *     | 
    | Divisão              |     /     |  
    | Soma                 |     +     | 
    | Subtração            |     -     |
    | Potenciação          |    **     |   
    | Raiz Quadrada        |   sqrt(x) | 
    | ln                   |   log(x)  |  
    | log base 10          |  log10(x) |
    | log base 2           |   log2(x) |
    | ln                   |   log(x)  |
    | Seno                 |   sin(x)  | 
    | Cosseno              |   cos(x)  |  
    | Tangente             |   tan(x)  | 
    | Seno hiperbólico     |   sinh(x) | 
    | Cosseno hiperbólico  |   cosh(x) |  
    | Tangente hiperbólica |   tanh(x) |   
    | Arco Seno            |   asin(x) | 
    | Arco Cosseno         |   acos(x) |  
    | Arco Tangente        |   atan(x) | 
    | e(Num de euler)      |      e    | 
    | pi                   |     pi    | 
    

    Caso a função desejada não possa ser feita com elementos da tabela, assuma que 
    o programa não será capaz de executa-la.

    Exemplos de funções corretamente inseridas:

        - (3*log(x)+4)
        - (e**(x**2))
  
    Exemplos de funções inseridas de maneira INCORRETA:

        - 3log(x)
        - (eˆ(xˆ2))
    


    OBSERVAÇÕES ADICIONAIS

        >> Separador decimal é o PONTO.
        >> Valor negativo basta colocar - na frente.
        >> EX: -9.23456


## Dados de saída
Depois que você inseriu seus dados, o programa vai retornar os valores solicitados por você, 
isto é, vai receber os vetores x e y, sendo x_n as soluções do sistema.