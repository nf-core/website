---
title: nf-core/modules new modules
category: components
slack: https://nfcore.slack.com/archives/C09LJTQQ3EY
location: UFRN
leaders:
  mribeirodantas:
    name: Marcel Ribeiro-Dantas
    slack: https://nfcore.slack.com/team/U03932BSX1V
---

Este projeto de hackathon tem como foco criar
[novos módulos no nf-core/modules](https://github.com/nf-core/modules),
seguindo os padrões oficiais da comunidade nf-core.

A criação de módulos amplia a cobertura de ferramentas reutilizáveis,
facilita a manutenção de pipelines e fortalece o ecossistema colaborativo do nf-core.

Recursos:

- [Guia oficial para criação de módulos](https://nf-co.re/docs/contributing/contribute-components)
- [Repositório nf-core/modules](https://github.com/nf-core/modules)
- [Hello nf-core](https://training.nextflow.io/latest/hello_nf-core/)

---

## Objetivo

Criar novos módulos no nf-core/modules com testes,
documentação e lint em conformidade com os padrões do projeto.

Este trabalho está organizado para que cada participante implemente
um módulo por vez, tornando a contribuição acessível para iniciantes.

---

## O que os participantes farão

Cada colaborador irá:

1. Escolher uma ferramenta que ainda não possua módulo no repositório.
2. Verificar se já existe issue aberta para essa ferramenta.
3. Criar o módulo seguindo o template e as diretrizes do nf-core.
4. Adicionar testes, metadados e documentação.
5. Executar lint e testes localmente.
6. Abrir um Pull Request para revisão.

Cada novo módulo deve estar vinculado a uma issue no repositório.

---

## Tarefas

### Criar um novo módulo no nf-core/modules

1. Escolha uma ferramenta de bioinformática que ainda não tenha módulo no [nf-core/modules](https://nf-co.re/modules).

2. Abra uma issue (ou assuma uma issue existente) para evitar trabalho duplicado.

3. Siga as diretrizes do [guia oficial de contribuição de módulos](https://nf-co.re/docs/contributing/contribute-components).

4. Implemente o módulo incluindo:
   - Script/module file com entradas e saídas padronizadas
   - Metadados e documentação de uso
   - Testes mínimos para validação funcional

5. Rode os comandos de lint e testes localmente.

6. Abra um Pull Request referenciando a issue correspondente.

7. Trate feedback de CI e revisão até a mesclagem.

Após concluir, os colaboradores são encorajados a criar módulos adicionais.

---

## Preparação recomendada

Idealmente, os participantes devem ter:

- Familiaridade básica com Git e GitHub (fork, branch e pull request)
- Conhecimento básico de Nextflow
- Familiaridade com estruturas do nf-core/modules

O seguinte material de treinamento é recomendado:

- [Hello nextflow](https://training.nextflow.io/latest/hello_nextflow/)
- [Hello nf-core](https://training.nextflow.io/latest/hello_nf-core/)
