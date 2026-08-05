---
title: nf-core/modules topics migration
category: components
slack: https://nfcore.slack.com/archives/C09LJTQQ3EY
location: UFRN
leaders:
  mribeirodantas:
    name: Marcel Ribeiro-Dantas
    slack: https://nfcore.slack.com/team/U03932BSX1V
---

Este projeto de hackathon tem como foco atualizar
[módulos do nf-core](https://github.com/nf-core/modules) para usar recursos modernos
do Nextflow, especificamente **topic channels**.

_Topic channels_ simplificam a coleta de versões entre módulos
e estão se tornando a abordagem padrão dentro do nf-core.

Recursos:

- [Tutorial de migração para topics do nf-core](https://nf-co.re/docs/tutorials/migrate_to_topics/update_modules)
- [Topic channels na documentação do Nextflow](https://nextflow.io/docs/latest/reference/channel.html#topic)
- [Topic channels no dashboard de estatísticas de adoção do nf-core](https://nf-core-stats.netlify.app/code/container_conversion/#version-topics-adoption-over-time)

---

## Objetivo

Migrar todos os módulos do nf-core para usar topic channels.

Este trabalho está organizado para que cada participante atualize um módulo por vez,
o que o torna uma contribuição amigável para iniciantes.

---

## O que os participantes farão

Cada colaborador irá:

1. Escolher um módulo que ainda não tenha sido migrado.
2. Atualizá-lo para usar topic channels seguindo as diretrizes oficiais.
3. Executar testes e verificações de lint localmente.
4. Abrir um Pull Request para revisão.

Cada migração corresponde a uma issue no repositório de módulos.

---

## Tarefas

### Migrar um módulo para topic channels

1. Escolha um módulo disponível no [rastreador de migração](https://github.com/nf-core/modules/issues/9978)

2. Atribua a issue a você para evitar trabalho duplicado.

3. Siga o [guia oficial de migração](https://nf-co.re/docs/tutorials/migrate_to_topics/update_modules). Estatísticas [aqui](https://nf-co.re/stats/code/container_conversion/).

4. Atualize o módulo:
  - Substituindo a lógica legada de coleta de versões
  - Usando topic channels para reporte de versões
  - Atualizando os testes, se necessário

5. Execute os testes do módulo e as verificações de lint localmente.

6. Abra um Pull Request referenciando a issue de migração.

7. Trate feedback de CI ou de revisão até a mesclagem.

Após concluir, os colaboradores são encorajados a migrar módulos adicionais.

---

## Preparação recomendada

Idealmente, os participantes devem ter:

- Familiaridade básica com Git e GitHub (fazer fork de repositórios, criar branches e abrir pull requests)
- Conhecimento básico de Nextflow
- Familiaridade com módulos nf-core

O seguinte material de treinamento é recomendado:

- [Hello nextflow](https://training.nextflow.io/latest/hello_nextflow/)
- [Hello nf-core](https://training.nextflow.io/latest/hello_nf-core/)
